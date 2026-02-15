#include "manual_stitch_page.hpp"

#include "gapneedle/fasta_io.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QCompleter>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QTextEdit>
#include <QVBoxLayout>

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace {

QComboBox* createSearchableCombo(QWidget* parent) {
  auto* combo = new QComboBox(parent);
  combo->setEditable(true);
  combo->setInsertPolicy(QComboBox::NoInsert);
  auto* c = combo->completer();
  c->setCompletionMode(QCompleter::PopupCompletion);
  c->setFilterMode(Qt::MatchContains);
  c->setCaseSensitivity(Qt::CaseInsensitive);
  return combo;
}

QString normalizedFsPath(QString path) {
  path = path.trimmed();
  if ((path.startsWith('"') && path.endsWith('"')) ||
      (path.startsWith('\'') && path.endsWith('\''))) {
    path = path.mid(1, path.size() - 2).trimmed();
  }
  if (path.isEmpty()) return {};
  return QDir::fromNativeSeparators(path);
}

QString sourceTitle(const QString& key) {
  if (key == "t") {
    return "Target";
  }
  if (key == "q") {
    return "Query";
  }
  return key;
}

}  // namespace

ManualStitchPage::ManualStitchPage(gapneedle::GapNeedleFacade* facade, QWidget* parent)
    : QWidget(parent), facade_(facade) {
  auto* root = new QVBoxLayout(this);
  root->setContentsMargins(12, 12, 12, 12);
  root->setSpacing(8);

  auto* header = new QHBoxLayout();
  auto* title = new QLabel("Manual stitch (coordinate mode)", this);
  QFont tf = title->font();
  tf.setBold(true);
  title->setFont(tf);
  header->addWidget(title);

  auto* loadLogBtn = new QPushButton("Load log (.json/.md)", this);
  connect(loadLogBtn, &QPushButton::clicked, this, &ManualStitchPage::onLoadLog);
  header->addWidget(loadLogBtn);
  header->addStretch(1);
  root->addLayout(header);

  auto* form = new QFormLayout();
  targetFasta_ = new QLineEdit(this);
  queryFasta_ = new QLineEdit(this);
  pafPath_ = new QLineEdit(this);
  pafPath_->setReadOnly(true);
  contextSpin_ = new QSpinBox(this);
  contextSpin_->setRange(20, 5000);
  contextSpin_->setValue(200);

  auto* targetRow = new QWidget(this);
  auto* targetLayout = new QHBoxLayout(targetRow);
  targetLayout->setContentsMargins(0, 0, 0, 0);
  auto* browseTarget = new QPushButton("Browse", targetRow);
  auto* loadTarget = new QPushButton("Load names", targetRow);
  targetLayout->addWidget(targetFasta_, 1);
  targetLayout->addWidget(browseTarget);
  targetLayout->addWidget(loadTarget);

  auto* queryRow = new QWidget(this);
  auto* queryLayout = new QHBoxLayout(queryRow);
  queryLayout->setContentsMargins(0, 0, 0, 0);
  auto* browseQuery = new QPushButton("Browse", queryRow);
  auto* loadQuery = new QPushButton("Load names", queryRow);
  queryLayout->addWidget(queryFasta_, 1);
  queryLayout->addWidget(browseQuery);
  queryLayout->addWidget(loadQuery);

  form->addRow("Target FASTA", targetRow);
  form->addRow("Query FASTA", queryRow);
  form->addRow("PAF (from Align)", pafPath_);
  form->addRow("Breakpoint context", contextSpin_);
  root->addLayout(form);

  connect(browseTarget, &QPushButton::clicked, this, &ManualStitchPage::onBrowseTarget);
  connect(browseQuery, &QPushButton::clicked, this, &ManualStitchPage::onBrowseQuery);
  connect(loadTarget, &QPushButton::clicked, this, &ManualStitchPage::onLoadTargetNames);
  connect(loadQuery, &QPushButton::clicked, this, &ManualStitchPage::onLoadQueryNames);

  auto* extraHeader = new QHBoxLayout();
  extraHeader->addWidget(new QLabel("Extra FASTA sources", this));
  auto* addExtraBtn = new QPushButton("Add FASTA", this);
  addExtraBtn->setObjectName("primaryButton");
  connect(addExtraBtn, &QPushButton::clicked, this, &ManualStitchPage::onAddExtraSource);
  extraHeader->addWidget(addExtraBtn);
  extraHeader->addStretch(1);
  root->addLayout(extraHeader);

  auto* extraContainer = new QWidget(this);
  extraRowsLayout_ = new QVBoxLayout(extraContainer);
  extraRowsLayout_->setContentsMargins(0, 0, 0, 0);
  extraRowsLayout_->setSpacing(4);
  root->addWidget(extraContainer);

  auto* segRow = new QHBoxLayout();
  sourceCombo_ = new QComboBox(this);
  sourceCombo_->setMinimumWidth(120);
  seqCombo_ = createSearchableCombo(this);
  startEdit_ = new QLineEdit(this);
  endEdit_ = new QLineEdit(this);
  startEdit_->setPlaceholderText("start");
  endEdit_->setPlaceholderText("end");
  auto* reverseCheck = new QCheckBox("Reverse-complement", this);
  reverseBtn_ = new QPushButton("Add segment", this);
  reverseBtn_->setObjectName("primaryButton");
  auto* resumeBtn = new QPushButton("Resume", this);

  segRow->addWidget(sourceCombo_);
  segRow->addWidget(seqCombo_, 2);
  segRow->addWidget(startEdit_);
  segRow->addWidget(endEdit_);
  segRow->addWidget(reverseCheck);
  segRow->addWidget(reverseBtn_);
  segRow->addWidget(resumeBtn);
  segRow->addStretch(1);
  root->addLayout(segRow);

  connect(sourceCombo_, &QComboBox::currentTextChanged, this, &ManualStitchPage::onSourceChanged);
  connect(reverseBtn_, &QPushButton::clicked, this, [this, reverseCheck]() {
    reverseBtn_->setProperty("reverse_checked", reverseCheck->isChecked());
    onAddSegment();
  });
  connect(resumeBtn, &QPushButton::clicked, this, &ManualStitchPage::onResumeSegment);

  segmentList_ = new QListWidget(this);
  connect(segmentList_, &QListWidget::itemSelectionChanged, this, &ManualStitchPage::onSegmentSelectionChanged);
  root->addWidget(segmentList_, 1);

  auto* detailRow = new QHBoxLayout();
  preview_ = new QTextEdit(this);
  preview_->setReadOnly(true);
  preview_->setMinimumHeight(220);
  detail_ = new QTextEdit(this);
  detail_->setReadOnly(true);
  detail_->setMinimumWidth(320);
  detailRow->addWidget(preview_, 2);
  detailRow->addWidget(detail_, 1);
  root->addLayout(detailRow);

  auto* footer = new QHBoxLayout();
  auto* checkBtn = new QPushButton("Check breakpoints", this);
  auto* exportBtn = new QPushButton("Export merged FASTA", this);
  exportBtn->setObjectName("primaryButton");
  auto* upBtn = new QPushButton("Move up", this);
  auto* downBtn = new QPushButton("Move down", this);
  auto* removeBtn = new QPushButton("Remove selected", this);
  footer->addWidget(checkBtn);
  footer->addWidget(exportBtn);
  footer->addWidget(upBtn);
  footer->addWidget(downBtn);
  footer->addWidget(removeBtn);
  footer->addStretch(1);
  root->addLayout(footer);

  connect(checkBtn, &QPushButton::clicked, this, &ManualStitchPage::onCheckBreakpoints);
  connect(exportBtn, &QPushButton::clicked, this, &ManualStitchPage::onExport);
  connect(upBtn, &QPushButton::clicked, this, &ManualStitchPage::onMoveSegmentUp);
  connect(downBtn, &QPushButton::clicked, this, &ManualStitchPage::onMoveSegmentDown);
  connect(removeBtn, &QPushButton::clicked, this, &ManualStitchPage::onRemoveSegment);

  result_ = new QTextEdit(this);
  result_->setReadOnly(true);
  result_->setMinimumHeight(120);
  root->addWidget(result_);

  namesBySource_["t"] = {};
  namesBySource_["q"] = {};
  refreshSourceCombo("t");
}

void ManualStitchPage::setAlignmentContext(const QString& targetFasta,
                                           const QString& queryFasta,
                                           const QString& targetSeq,
                                           const QString& querySeq,
                                           const QString& pafPath) {
  targetFasta_->setText(targetFasta);
  queryFasta_->setText(queryFasta);
  pafPath_->setText(pafPath);

  loadNamesForSource("t", targetFasta);
  loadNamesForSource("q", queryFasta);

  if (!targetSeq.isEmpty()) {
    seqCombo_->setCurrentText(targetSeq);
  }

  appendResult(QString("Alignment context synced: target=%1 query=%2").arg(targetSeq, querySeq));
}

void ManualStitchPage::onBrowseTarget() {
  const QString path = QFileDialog::getOpenFileName(this, "Select target FASTA", QString(),
                                                    "FASTA (*.fa *.fasta *.fna);;All files (*)");
  if (path.isEmpty()) {
    return;
  }
  targetFasta_->setText(path);
  loadNamesForSource("t", path);
  if (sourceCombo_->currentData().toString() == "t") {
    refreshSeqCombo();
  }
}

void ManualStitchPage::onBrowseQuery() {
  const QString path = QFileDialog::getOpenFileName(this, "Select query FASTA", QString(),
                                                    "FASTA (*.fa *.fasta *.fna);;All files (*)");
  if (path.isEmpty()) {
    return;
  }
  queryFasta_->setText(path);
  loadNamesForSource("q", path);
  if (sourceCombo_->currentData().toString() == "q") {
    refreshSeqCombo();
  }
}

void ManualStitchPage::onLoadTargetNames() {
  if (loadNamesForSource("t", targetFasta_->text().trimmed(), true) &&
      sourceCombo_->currentData().toString() == "t") {
    refreshSeqCombo();
  }
}

void ManualStitchPage::onLoadQueryNames() {
  if (loadNamesForSource("q", queryFasta_->text().trimmed(), true) &&
      sourceCombo_->currentData().toString() == "q") {
    refreshSeqCombo();
  }
}

void ManualStitchPage::onAddExtraSource() {
  const QString key = QString("x%1").arg(nextExtraId_++);
  auto* row = new QWidget(this);
  auto* rowLayout = new QHBoxLayout(row);
  rowLayout->setContentsMargins(0, 0, 0, 0);

  auto* label = new QLabel(QString("%1 FASTA").arg(key), row);
  auto* path = new QLineEdit(row);
  auto* browse = new QPushButton("Browse", row);
  auto* load = new QPushButton("Load names", row);
  auto* remove = new QPushButton("Remove", row);

  rowLayout->addWidget(label);
  rowLayout->addWidget(path, 1);
  rowLayout->addWidget(browse);
  rowLayout->addWidget(load);
  rowLayout->addWidget(remove);
  extraRowsLayout_->addWidget(row);

  ExtraSourceRow ext;
  ext.key = key;
  ext.rowWidget = row;
  ext.pathEdit = path;
  extras_[key] = ext;
  namesBySource_[key] = {};

  connect(browse, &QPushButton::clicked, this, [this, key, path]() {
    const QString p = QFileDialog::getOpenFileName(this,
                                                   QString("Select FASTA for %1").arg(key),
                                                   QString(),
                                                   "FASTA (*.fa *.fasta *.fna);;All files (*)");
    if (!p.isEmpty()) {
      path->setText(p);
    }
  });
  connect(load, &QPushButton::clicked, this, [this, key, path]() {
    if (loadNamesForSource(key, path->text().trimmed(), true) &&
        sourceCombo_->currentData().toString() == key) {
      refreshSeqCombo();
    }
  });
  connect(remove, &QPushButton::clicked, this, [this, key]() {
    for (const auto& seg : segments_) {
      if (seg.source == key) {
        QMessageBox::warning(this, "Source in use", "Remove segments using this source first.");
        return;
      }
    }
    auto it = extras_.find(key);
    if (it != extras_.end()) {
      it->rowWidget->deleteLater();
      extras_.erase(it);
    }
    namesBySource_.remove(key);
    refreshSourceCombo("t");
    refreshSeqCombo();
  });

  refreshSourceCombo(key);
}

void ManualStitchPage::onSourceChanged() { refreshSeqCombo(); }

void ManualStitchPage::onAddSegment() {
  const QString source = sourceCombo_->currentData().toString();
  const QString name = seqCombo_->currentText().trimmed();
  bool okStart = false;
  bool okEnd = false;
  const int start = startEdit_->text().replace(",", "").toInt(&okStart);
  const int end = endEdit_->text().replace(",", "").toInt(&okEnd);
  const bool reverse = reverseBtn_->property("reverse_checked").toBool();

  if (source.isEmpty() || name.isEmpty() || !okStart || !okEnd || end <= start) {
    QMessageBox::warning(this, "Invalid segment", "Please set source, sequence, and valid start/end.");
    return;
  }

  SegmentItem item;
  item.source = source;
  item.seqName = name;
  item.start = start;
  item.end = end;
  item.reverse = reverse;
  segments_.push_back(item);

  startEdit_->clear();
  endEdit_->clear();
  refreshSegments();
}

void ManualStitchPage::onRemoveSegment() {
  const int row = segmentList_->currentRow();
  if (row < 0 || row >= static_cast<int>(segments_.size())) {
    return;
  }
  segments_.erase(segments_.begin() + row);
  refreshSegments();
}

void ManualStitchPage::onMoveSegmentUp() {
  const int row = segmentList_->currentRow();
  if (row <= 0 || row >= static_cast<int>(segments_.size())) {
    return;
  }
  std::swap(segments_[row], segments_[row - 1]);
  refreshSegments();
  segmentList_->setCurrentRow(row - 1);
}

void ManualStitchPage::onMoveSegmentDown() {
  const int row = segmentList_->currentRow();
  if (row < 0 || row + 1 >= static_cast<int>(segments_.size())) {
    return;
  }
  std::swap(segments_[row], segments_[row + 1]);
  refreshSegments();
  segmentList_->setCurrentRow(row + 1);
}

void ManualStitchPage::onResumeSegment() {
  const int row = segmentList_->currentRow();
  if (row < 0 || row >= static_cast<int>(segments_.size())) {
    return;
  }
  const auto& seg = segments_[row];
  const int idx = sourceCombo_->findData(seg.source);
  if (idx >= 0) {
    sourceCombo_->setCurrentIndex(idx);
  }
  seqCombo_->setCurrentText(seg.seqName);
  startEdit_->setText(QString::number(seg.start));
  endEdit_->setText(QString::number(seg.end));
  reverseBtn_->setProperty("reverse_checked", seg.reverse);
}

void ManualStitchPage::onSegmentSelectionChanged() {
  const int row = segmentList_->currentRow();
  if (row < 0 || row >= static_cast<int>(segments_.size())) {
    detail_->clear();
    return;
  }
  const auto& seg = segments_[row];
  QString text = QString("Segment [%1] %2:%3 %4-%5%6\n")
                     .arg(row)
                     .arg(seg.source)
                     .arg(seg.seqName)
                     .arg(seg.start)
                     .arg(seg.end)
                     .arg(seg.reverse ? " (RC)" : "");
  if (!seg.seq.isEmpty()) {
    text += "\nLeft context:\n" + seg.leftBefore + "\n" + seg.leftAfter + "\n";
    text += "\nRight context:\n" + seg.rightBefore + "\n" + seg.rightAfter + "\n";
  } else {
    text += "\nSegment not materialized yet. Click 'Check breakpoints' or export.\n";
  }
  detail_->setPlainText(text);
}

void ManualStitchPage::onCheckBreakpoints() {
  if (segments_.empty()) {
    QMessageBox::information(this, "No segments", "Add at least one segment.");
    return;
  }
  if (!materializeAll(contextSpin_->value())) {
    return;
  }
  refreshPreview();
  appendResult(allBreakpointsMatch() ? "All breakpoints passed." : "Breakpoint differences detected.");
}

void ManualStitchPage::onExport() {
  if (segments_.empty()) {
    QMessageBox::information(this, "No segments", "Add at least one segment.");
    return;
  }
  if (!materializeAll(contextSpin_->value())) {
    return;
  }

  QString out = QFileDialog::getSaveFileName(this,
                                             "Save merged FASTA",
                                             QFileInfo(targetFasta_->text()).absolutePath() + "/stitched.fa",
                                             "FASTA (*.fa *.fasta)");
  if (out.isEmpty()) {
    return;
  }

  QString merged;
  for (const auto& seg : segments_) {
    merged += seg.seq;
  }

  gapneedle::FastaMap records;
  records[QFileInfo(out).completeBaseName().toStdString()] = merged.toStdString();
  gapneedle::writeFasta(out.toStdString(), records);

  QJsonObject root;
  root["target_fasta"] = targetFasta_->text();
  root["query_fasta"] = queryFasta_->text();
  root["paf"] = pafPath_->text();
  root["context_bp"] = contextSpin_->value();
  QJsonArray segs;
  for (const auto& seg : segments_) {
    QJsonObject s;
    s["source"] = seg.source;
    s["name"] = seg.seqName;
    s["start"] = seg.start;
    s["end"] = seg.end;
    s["reverse"] = seg.reverse;
    segs.append(s);
  }
  root["segments"] = segs;

  const QString logPath = out + ".session.json";
  QFile f(logPath);
  if (f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    f.close();
  }

  appendResult(QString("Saved: %1").arg(out));
  appendResult(QString("Session: %1").arg(logPath));
}

void ManualStitchPage::onLoadLog() {
  const QString path = QFileDialog::getOpenFileName(this,
                                                    "Load stitch log",
                                                    QString(),
                                                    "Session (*.json *.md);;All files (*)");
  if (path.isEmpty()) {
    return;
  }

  QFile f(path);
  if (!f.open(QIODevice::ReadOnly)) {
    QMessageBox::critical(this, "Read failed", "Cannot read selected log file.");
    return;
  }
  const QByteArray raw = f.readAll();
  f.close();

  QJsonParseError err;
  const QJsonDocument doc = QJsonDocument::fromJson(raw, &err);
  if (err.error == QJsonParseError::NoError && doc.isObject()) {
    const QJsonObject root = doc.object();
    targetFasta_->setText(root.value("target_fasta").toString());
    queryFasta_->setText(root.value("query_fasta").toString());
    pafPath_->setText(root.value("paf").toString());
    contextSpin_->setValue(root.value("context_bp").toInt(200));

    loadNamesForSource("t", targetFasta_->text());
    loadNamesForSource("q", queryFasta_->text());

    segments_.clear();
    for (const auto& v : root.value("segments").toArray()) {
      const QJsonObject s = v.toObject();
      SegmentItem seg;
      seg.source = s.value("source").toString();
      seg.seqName = s.value("name").toString();
      seg.start = s.value("start").toInt();
      seg.end = s.value("end").toInt();
      seg.reverse = s.value("reverse").toBool(false);
      segments_.push_back(seg);
    }

    refreshSourceCombo("t");
    refreshSegments();
    appendResult(QString("Session loaded: %1").arg(path));
    return;
  }

  // Fallback parser for legacy markdown session logs.
  const QString text = QString::fromUtf8(raw);
  QRegularExpression lineRe(
      R"(^-\s*\[\d+\]\s*[^()]*\((t|q|x\d+)\)\s*(.+?)\s+(\d+)\s*-\s*(\d+)\s+.*$)");
  QRegularExpression targetRe(R"(^-\s*(Target FASTA|目标 FASTA):\s*(.+)$)");
  QRegularExpression queryRe(R"(^-\s*(Query FASTA|查询 FASTA):\s*(.+)$)");
  QRegularExpression ctxRe(R"((each|各取)\s*(\d+)\s*bp)");

  segments_.clear();
  for (const QString& rawLine : text.split('\n')) {
    const QString line = rawLine.trimmed();
    auto tm = targetRe.match(line);
    if (tm.hasMatch()) {
      targetFasta_->setText(tm.captured(2).trimmed());
      continue;
    }
    auto qm = queryRe.match(line);
    if (qm.hasMatch()) {
      queryFasta_->setText(qm.captured(2).trimmed());
      continue;
    }
    auto cm = ctxRe.match(line);
    if (cm.hasMatch()) {
      contextSpin_->setValue(cm.captured(2).toInt());
      continue;
    }

    auto m = lineRe.match(line);
    if (!m.hasMatch()) {
      continue;
    }
    SegmentItem seg;
    seg.source = m.captured(1).trimmed();
    QString name = m.captured(2).trimmed();
    if (name.endsWith("(RC)")) {
      name.chop(4);
      name = name.trimmed();
      seg.reverse = true;
    }
    seg.seqName = name;
    seg.start = m.captured(3).toInt();
    seg.end = m.captured(4).toInt();
    segments_.push_back(seg);
  }

  if (segments_.empty()) {
    QMessageBox::warning(this, "Unsupported log", "Cannot parse session from this file.");
    return;
  }
  loadNamesForSource("t", targetFasta_->text());
  loadNamesForSource("q", queryFasta_->text());
  refreshSourceCombo("t");
  refreshSegments();
  appendResult(QString("Legacy markdown session loaded: %1").arg(path));
}

void ManualStitchPage::appendResult(const QString& text) { result_->append(text); }

bool ManualStitchPage::loadNamesForSource(const QString& sourceKey, const QString& fastaPath, bool verbose) {
  const QString path = normalizedFsPath(fastaPath);
  if (path.isEmpty()) {
    namesBySource_[sourceKey] = {};
    if (verbose) {
      appendResult(QString("%1: FASTA path is empty.").arg(sourceTitle(sourceKey)));
    }
    return false;
  }

  const QFileInfo fi(path);
  if (!fi.exists() || !fi.isFile()) {
    namesBySource_[sourceKey] = {};
    if (verbose) {
      appendResult(QString("%1: FASTA not found: %2").arg(sourceTitle(sourceKey), path));
      QMessageBox::warning(this, "Load names failed",
                           QString("Cannot find FASTA file:\n%1").arg(path));
    }
    return false;
  }

  const QStringList names = fastaNamesFast(path);
  namesBySource_[sourceKey] = names;
  if (verbose) {
    if (names.isEmpty()) {
      appendResult(QString("%1: no sequence headers found in %2").arg(sourceTitle(sourceKey), path));
      QMessageBox::warning(this, "Load names failed",
                           QString("No FASTA headers were parsed from:\n%1").arg(path));
      return false;
    }
    appendResult(QString("%1: loaded %2 sequence names from %3")
                     .arg(sourceTitle(sourceKey))
                     .arg(names.size())
                     .arg(path));
  }
  return !names.isEmpty();
}

QString ManualStitchPage::sourcePath(const QString& sourceKey) const {
  if (sourceKey == "t") {
    return normalizedFsPath(targetFasta_->text());
  }
  if (sourceKey == "q") {
    return normalizedFsPath(queryFasta_->text());
  }
  auto it = extras_.find(sourceKey);
  if (it != extras_.end() && it->pathEdit) {
    return normalizedFsPath(it->pathEdit->text());
  }
  return {};
}

QStringList ManualStitchPage::sourceNames(const QString& sourceKey) const {
  auto it = namesBySource_.find(sourceKey);
  if (it == namesBySource_.end()) {
    return {};
  }
  return it.value();
}

void ManualStitchPage::refreshSourceCombo(const QString& keepKey) {
  QString current = keepKey;
  if (current.isEmpty()) {
    current = sourceCombo_->currentData().toString();
  }

  sourceCombo_->blockSignals(true);
  sourceCombo_->clear();
  sourceCombo_->addItem("t (target)", "t");
  sourceCombo_->addItem("q (query)", "q");
  for (auto it = extras_.cbegin(); it != extras_.cend(); ++it) {
    sourceCombo_->addItem(QString("%1 (extra)").arg(it.key()), it.key());
  }
  int idx = sourceCombo_->findData(current);
  if (idx < 0) {
    idx = 0;
  }
  sourceCombo_->setCurrentIndex(idx);
  sourceCombo_->blockSignals(false);
}

void ManualStitchPage::refreshSeqCombo() {
  const QString src = sourceCombo_->currentData().toString();
  const QString current = seqCombo_->currentText();
  const QStringList names = sourceNames(src);
  seqCombo_->blockSignals(true);
  seqCombo_->clear();
  seqCombo_->addItems(names);
  if (!current.isEmpty() && names.contains(current)) {
    seqCombo_->setCurrentText(current);
  }
  seqCombo_->blockSignals(false);
}

void ManualStitchPage::refreshSegments() {
  segmentList_->clear();
  long long total = 0;
  for (int i = 0; i < static_cast<int>(segments_.size()); ++i) {
    const auto& seg = segments_[i];
    const int len = seg.end - seg.start;
    total += len;
    segmentList_->addItem(QString("[%1] %2:%3 %4-%5 %6bp%7")
                              .arg(i)
                              .arg(seg.source)
                              .arg(seg.seqName)
                              .arg(seg.start)
                              .arg(seg.end)
                              .arg(len)
                              .arg(seg.reverse ? " (RC)" : ""));
  }
  appendResult(QString("Segments: %1 total length: %2").arg(segments_.size()).arg(total));
  refreshPreview();
}

void ManualStitchPage::refreshPreview() {
  preview_->clear();
  if (segments_.size() < 2) {
    preview_->setPlainText("Add at least two segments to preview breakpoints.");
    return;
  }
  bool hasContext = true;
  for (const auto& seg : segments_) {
    if (seg.rightBefore.isEmpty() && seg.leftAfter.isEmpty()) {
      hasContext = false;
      break;
    }
  }
  if (!hasContext) {
    preview_->setPlainText("Run 'Check breakpoints' to materialize context.");
    return;
  }

  for (int i = 0; i + 1 < static_cast<int>(segments_.size()); ++i) {
    const auto& left = segments_[i];
    const auto& right = segments_[i + 1];
    preview_->append(QString("[%1] %2:%3 %4-%5 -> %6:%7 %8-%9")
                         .arg(i)
                         .arg(left.source)
                         .arg(left.seqName)
                         .arg(left.start)
                         .arg(left.end)
                         .arg(right.source)
                         .arg(right.seqName)
                         .arg(right.start)
                         .arg(right.end));
    preview_->append(junctionPreview(left.rightBefore, right.leftAfter, contextSpin_->value()));

    const int ll = std::min(50, std::min(static_cast<int>(left.rightBefore.size()),
                                         static_cast<int>(right.leftBefore.size())));
    if (ll > 0) {
      const bool ok = left.rightBefore.right(ll) == right.leftBefore.right(ll);
      preview_->append(ok ? "Left flanks match" : "Left flanks differ");
    }

    const int rl = std::min(50, std::min(static_cast<int>(left.rightAfter.size()),
                                         static_cast<int>(right.leftAfter.size())));
    if (rl > 0) {
      const bool ok = left.rightAfter.left(rl) == right.leftAfter.left(rl);
      preview_->append(ok ? "Right flanks match" : "Right flanks differ");
    }
    preview_->append("");
  }
}

bool ManualStitchPage::materializeAll(int contextBp) {
  for (auto& seg : segments_) {
    if (!materializeSegment(seg, contextBp)) {
      return false;
    }
  }
  return true;
}

bool ManualStitchPage::materializeSegment(SegmentItem& seg, int contextBp) {
  try {
    seg.seq = readSegment(seg.source, seg.seqName, seg.start, seg.end, seg.reverse);
    seg.leftBefore = readSegment(seg.source, seg.seqName, std::max(0, seg.start - contextBp), seg.start, seg.reverse);
    seg.leftAfter = readSegment(seg.source, seg.seqName, seg.start, seg.start + contextBp, seg.reverse);
    seg.rightBefore = readSegment(seg.source, seg.seqName, std::max(0, seg.end - contextBp), seg.end, seg.reverse);
    seg.rightAfter = readSegment(seg.source, seg.seqName, seg.end, seg.end + contextBp, seg.reverse);
    return true;
  } catch (const std::exception& e) {
    QMessageBox::critical(this, "Materialize failed", e.what());
    return false;
  }
}

QStringList ManualStitchPage::fastaNamesFast(const QString& path) const {
  QStringList names;
  std::ifstream in(path.toStdString());
  if (!in) {
    return names;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      QString name = QString::fromStdString(line.substr(1)).trimmed();
      int sp = name.indexOf(' ');
      if (sp > 0) {
        name = name.left(sp);
      }
      if (!name.isEmpty()) {
        names.push_back(name);
      }
    }
  }
  names.removeDuplicates();
  return names;
}

QString ManualStitchPage::readSegment(const QString& sourceKey,
                                      const QString& seqName,
                                      int start,
                                      int end,
                                      bool reverse) const {
  const QString path = sourcePath(sourceKey);
  if (path.isEmpty()) {
    throw std::runtime_error(("Missing FASTA path for source: " + sourceKey).toStdString());
  }
  const QFileInfo fi(path);
  if (!fi.exists() || !fi.isFile()) {
    throw std::runtime_error(("FASTA file not found for source " + sourceKey + ": " + path).toStdString());
  }
  auto it = fastaCache_.find(path.toStdString());
  if (it == fastaCache_.end()) {
    try {
      fastaCache_[path.toStdString()] = gapneedle::readFasta(path.toStdString());
    } catch (const std::exception& e) {
      throw std::runtime_error(("Failed to open FASTA for source " + sourceKey + ": " + path +
                                " (" + e.what() + ")").toStdString());
    }
    it = fastaCache_.find(path.toStdString());
  }

  auto sit = it->second.find(seqName.toStdString());
  if (sit == it->second.end()) {
    throw std::runtime_error(("Sequence not found: " + seqName).toStdString());
  }
  QString seq = QString::fromStdString(sit->second);
  if (reverse) {
    seq = QString::fromStdString(gapneedle::reverseComplement(seq.toStdString()));
  }
  start = std::max(0, start);
  end = std::min(end, static_cast<int>(seq.size()));
  if (end <= start) {
    return {};
  }
  return seq.mid(start, end - start);
}

QString ManualStitchPage::junctionPreview(const QString& left, const QString& right, int contextBp) const {
  const int ctx = std::max(0, contextBp);
  const QString l = left.right(std::min(ctx, static_cast<int>(left.size())));
  const QString r = right.left(std::min(ctx, static_cast<int>(right.size())));
  return l + "|" + r;
}

bool ManualStitchPage::allBreakpointsMatch() const {
  if (segments_.size() < 2) {
    return false;
  }
  for (int i = 0; i + 1 < static_cast<int>(segments_.size()); ++i) {
    const auto& a = segments_[i];
    const auto& b = segments_[i + 1];
    const int ll = std::min(50, std::min(static_cast<int>(a.rightBefore.size()),
                                         static_cast<int>(b.leftBefore.size())));
    const int rl = std::min(50, std::min(static_cast<int>(a.rightAfter.size()),
                                         static_cast<int>(b.leftAfter.size())));
    if (ll == 0 || rl == 0) {
      return false;
    }
    if (a.rightBefore.right(ll) != b.leftBefore.right(ll)) {
      return false;
    }
    if (a.rightAfter.left(rl) != b.leftAfter.left(rl)) {
      return false;
    }
  }
  return true;
}
