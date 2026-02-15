#include "align_page.hpp"

#include "gapneedle/fasta_io.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QCompleter>
#include <QCoreApplication>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFutureWatcher>
#include <QFormLayout>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QTextCharFormat>
#include <QTextCursor>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QtConcurrent/QtConcurrent>

#include <filesystem>

namespace {

QString safePart(QString in) {
  static const QRegularExpression bad(R"([^A-Za-z0-9._-])");
  in.replace(bad, "_");
  return in;
}

QComboBox* createSearchableCombo(QWidget* parent) {
  auto* combo = new QComboBox(parent);
  combo->setEditable(true);
  combo->setInsertPolicy(QComboBox::NoInsert);
  combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  combo->setMinimumContentsLength(24);
  auto* completer = combo->completer();
  completer->setFilterMode(Qt::MatchContains);
  completer->setCaseSensitivity(Qt::CaseInsensitive);
  return combo;
}

QString logLevelForKey(const QString& key) {
  if (key == "Config") return "config";
  if (key == "Cache hit" || key == "Cache miss") return "cache";
  if (key == "PAF" || key == "Status") return "success";
  if (key == "Error") return "error";
  if (key == "Running" || key == "Load") return "running";
  return "info";
}

QTextCharFormat badgeFormatForLevel(const QString& level) {
  QColor bg("#E9EAEE");
  QColor fg("#4A4A4A");
  if (level == "config") {
    bg = QColor("#E8F0FE");
    fg = QColor("#2557D6");
  } else if (level == "cache") {
    bg = QColor("#FFF4E5");
    fg = QColor("#A15A00");
  } else if (level == "running") {
    bg = QColor("#EAF2FF");
    fg = QColor("#3367D6");
  } else if (level == "success") {
    bg = QColor("#E8F7ED");
    fg = QColor("#1E7A3D");
  } else if (level == "error") {
    bg = QColor("#FDEBEC");
    fg = QColor("#B3261E");
  }
  QTextCharFormat fmt;
  fmt.setForeground(fg);
  fmt.setBackground(bg);
  fmt.setFontWeight(QFont::DemiBold);
  return fmt;
}

bool splitLogLine(const QString& text, QString* key, QString* body) {
  const int p = text.indexOf(':');
  if (p <= 0) {
    return false;
  }
  *key = text.left(p).trimmed();
  *body = text.mid(p + 1).trimmed();
  return !key->isEmpty();
}

}  // namespace

struct AlignTaskResult {
  gapneedle::AlignmentResult result;
  QString configLine;
  QString cacheLine;
  QString error;
};

struct FastaNameLoadResult {
  QStringList names;
  QString error;
};

AlignPage::AlignPage(gapneedle::GapNeedleFacade* facade, QWidget* parent)
    : QWidget(parent), facade_(facade) {
  auto* outer = new QVBoxLayout(this);
  outer->setContentsMargins(14, 14, 14, 14);
  outer->setSpacing(10);

  auto* card = new QFrame(this);
  card->setFrameShape(QFrame::StyledPanel);
  auto* layout = new QVBoxLayout(card);
  layout->setContentsMargins(14, 14, 14, 14);
  layout->setSpacing(10);

  auto* title = new QLabel("Alignment", card);
  QFont tf = title->font();
  tf.setPointSize(tf.pointSize() + 1);
  tf.setBold(true);
  title->setFont(tf);
  layout->addWidget(title);

  auto* form = new QFormLayout();
  form->setHorizontalSpacing(10);
  form->setVerticalSpacing(8);

  targetFasta_ = new QLineEdit(card);
  auto* targetRow = new QWidget(card);
  auto* targetRowLayout = new QHBoxLayout(targetRow);
  targetRowLayout->setContentsMargins(0, 0, 0, 0);
  targetRowLayout->setSpacing(6);
  auto* browseTarget = new QPushButton("Browse", targetRow);
  targetRowLayout->addWidget(targetFasta_, 1);
  targetRowLayout->addWidget(browseTarget);
  form->addRow("Target FASTA", targetRow);

  queryFasta_ = new QLineEdit(card);
  auto* queryRow = new QWidget(card);
  auto* queryRowLayout = new QHBoxLayout(queryRow);
  queryRowLayout->setContentsMargins(0, 0, 0, 0);
  queryRowLayout->setSpacing(6);
  auto* browseQuery = new QPushButton("Browse", queryRow);
  queryRowLayout->addWidget(queryFasta_, 1);
  queryRowLayout->addWidget(browseQuery);
  form->addRow("Query FASTA", queryRow);

  targetSeqCombo_ = createSearchableCombo(card);
  querySeqCombo_ = createSearchableCombo(card);
  form->addRow("Target sequence", targetSeqCombo_);
  form->addRow("Query sequence", querySeqCombo_);

  auto* optRow = new QWidget(card);
  auto* optLayout = new QHBoxLayout(optRow);
  optLayout->setContentsMargins(0, 0, 0, 0);
  optLayout->setSpacing(10);

  presetCombo_ = new QComboBox(optRow);
  presetCombo_->addItems({"asm5", "asm10", "asm20", "map-ont", "map-pb", "sr", "splice", "ava-ont", "ava-pb"});
  presetCombo_->setCurrentText("asm20");
  threads_ = new QSpinBox(optRow);
  threads_->setRange(1, 128);
  threads_->setValue(4);
  reverseTarget_ = new QCheckBox("Reverse-complement target", optRow);
  reverseQuery_ = new QCheckBox("Reverse-complement query", optRow);

  optLayout->addWidget(new QLabel("Preset", optRow));
  optLayout->addWidget(presetCombo_);
  optLayout->addWidget(new QLabel("Threads", optRow));
  optLayout->addWidget(threads_);
  optLayout->addWidget(reverseTarget_);
  optLayout->addWidget(reverseQuery_);
  optLayout->addStretch(1);
  form->addRow("Options", optRow);

  cachePathView_ = new QLineEdit(card);
  cachePathView_->setReadOnly(true);
  cachePathView_->setPlaceholderText("PAF cache path will be generated automatically");
  form->addRow("Cache PAF", cachePathView_);

  layout->addLayout(form);

  auto* btnRow = new QHBoxLayout();
  runBtn_ = new QPushButton("Run alignment", card);
  runBtn_->setObjectName("primaryButton");
  auto* clearBtn = new QPushButton("Clear log", card);
  btnRow->addWidget(runBtn_);
  btnRow->addWidget(clearBtn);
  btnRow->addStretch(1);
  layout->addLayout(btnRow);

  log_ = new QTextEdit(card);
  log_->setObjectName("alignLog");
  log_->setReadOnly(true);
  log_->setMinimumHeight(260);
  layout->addWidget(log_, 1);

  outer->addWidget(card, 1);

  connect(browseTarget, &QPushButton::clicked, this, &AlignPage::onBrowseTarget);
  connect(browseQuery, &QPushButton::clicked, this, &AlignPage::onBrowseQuery);
  connect(targetFasta_, &QLineEdit::editingFinished, this, &AlignPage::onTargetPathEdited);
  connect(queryFasta_, &QLineEdit::editingFinished, this, &AlignPage::onQueryPathEdited);
  connect(targetSeqCombo_, &QComboBox::currentTextChanged, this, [this]() { cachePathView_->setText(computeCachePafPath()); });
  connect(querySeqCombo_, &QComboBox::currentTextChanged, this, [this]() { cachePathView_->setText(computeCachePafPath()); });
  connect(presetCombo_, &QComboBox::currentTextChanged, this, [this]() { cachePathView_->setText(computeCachePafPath()); });
  connect(reverseTarget_, &QCheckBox::toggled, this, [this]() { cachePathView_->setText(computeCachePafPath()); });
  connect(reverseQuery_, &QCheckBox::toggled, this, [this]() { cachePathView_->setText(computeCachePafPath()); });
  connect(runBtn_, &QPushButton::clicked, this, &AlignPage::onRunAlign);
  connect(clearBtn, &QPushButton::clicked, log_, &QTextEdit::clear);
}

void AlignPage::appendLog(const QString& text) {
  QTextCursor c = log_->textCursor();
  c.movePosition(QTextCursor::End);

  QString key;
  QString body;
  if (!splitLogLine(text, &key, &body)) {
    c.insertText(text + "\n");
    log_->setTextCursor(c);
    return;
  }

  const QString level = logLevelForKey(key);
  QTextCharFormat badgeFmt = badgeFormatForLevel(level);
  badgeFmt.setFontPointSize(10.0);
  c.insertText(" " + key + " ", badgeFmt);

  QTextCharFormat bodyFmt;
  bodyFmt.setForeground(QColor("#3A3A3C"));
  bodyFmt.setFontPointSize(10.5);
  c.insertText(" " + body + "\n", bodyFmt);
  log_->setTextCursor(c);
}

void AlignPage::onBrowseTarget() {
  const QString path = QFileDialog::getOpenFileName(this,
                                                    "Select target FASTA",
                                                    QString(),
                                                    "FASTA (*.fa *.fasta *.fna);;All files (*)");
  if (path.isEmpty()) {
    return;
  }
  targetFasta_->setText(path);
  onTargetPathEdited();
}

void AlignPage::onBrowseQuery() {
  const QString path = QFileDialog::getOpenFileName(this,
                                                    "Select query FASTA",
                                                    QString(),
                                                    "FASTA (*.fa *.fasta *.fna);;All files (*)");
  if (path.isEmpty()) {
    return;
  }
  queryFasta_->setText(path);
  onQueryPathEdited();
}

void AlignPage::loadSequenceNames(const QString& fastaPath, QComboBox* combo) {
  const bool isTarget = combo == targetSeqCombo_;
  int* tokenPtr = isTarget ? &targetLoadToken_ : &queryLoadToken_;
  ++(*tokenPtr);
  const int token = *tokenPtr;

  combo->clear();
  if (fastaPath.trimmed().isEmpty()) {
    combo->setEnabled(true);
    return;
  }

  const QString key = normalizedPathKey(fastaPath);
  auto it = fastaNamesCache_.find(key);
  if (it != fastaNamesCache_.end()) {
    combo->addItems(it.value());
    combo->setEnabled(true);
    return;
  }

  combo->setEnabled(false);
  combo->addItem("Loading sequence names...");

  auto* watcher = new QFutureWatcher<FastaNameLoadResult>(this);
  connect(watcher, &QFutureWatcher<FastaNameLoadResult>::finished, this, [this, watcher, combo, key, isTarget, token]() {
    watcher->deleteLater();
    const int currentToken = isTarget ? targetLoadToken_ : queryLoadToken_;
    if (token != currentToken) {
      return;
    }

    FastaNameLoadResult loaded = watcher->result();
    combo->blockSignals(true);
    combo->clear();
    if (!loaded.error.isEmpty()) {
      combo->setEnabled(true);
      combo->blockSignals(false);
      appendLog(QString("Load: failed to load FASTA names: %1").arg(loaded.error));
      return;
    }

    fastaNamesCache_.insert(key, loaded.names);
    combo->addItems(loaded.names);
    combo->setEnabled(true);
    combo->blockSignals(false);
    cachePathView_->setText(computeCachePafPath());
  });

  const std::string path = fastaPath.toStdString();
  watcher->setFuture(QtConcurrent::run([path]() {
    FastaNameLoadResult out;
    try {
      const auto names = gapneedle::readFastaNames(path);
      out.names.reserve(static_cast<qsizetype>(names.size()));
      for (const auto& n : names) {
        out.names.push_back(QString::fromStdString(n));
      }
    } catch (const std::exception& e) {
      out.error = QString::fromUtf8(e.what());
    } catch (...) {
      out.error = "unknown FASTA read error";
    }
    return out;
  }));
}

QString AlignPage::normalizedPathKey(const QString& path) const {
  const QFileInfo fi(path.trimmed());
  const QString canonical = fi.canonicalFilePath();
  if (!canonical.isEmpty()) {
    return canonical;
  }
  const QString absolute = fi.absoluteFilePath();
  if (!absolute.isEmpty()) {
    return absolute;
  }
  return path.trimmed();
}

void AlignPage::onTargetPathEdited() {
  loadSequenceNames(targetFasta_->text().trimmed(), targetSeqCombo_);
  cachePathView_->setText(computeCachePafPath());
}

void AlignPage::onQueryPathEdited() {
  loadSequenceNames(queryFasta_->text().trimmed(), querySeqCombo_);
  cachePathView_->setText(computeCachePafPath());
}

QString AlignPage::computeCachePafPath() const {
  const QString tf = targetFasta_->text().trimmed();
  const QString qf = queryFasta_->text().trimmed();
  const QString ts = targetSeqCombo_->currentText().trimmed();
  const QString qs = querySeqCombo_->currentText().trimmed();
  if (tf.isEmpty() || qf.isEmpty() || ts.isEmpty() || qs.isEmpty()) {
    return {};
  }

  const QFileInfo tfi(tf);
  const QFileInfo qfi(qf);
  const QString tPart = safePart(tfi.completeBaseName()) + "." + safePart(ts) + (reverseTarget_->isChecked() ? "_rc" : "");
  const QString qPart = safePart(qfi.completeBaseName()) + "." + safePart(qs) + (reverseQuery_->isChecked() ? "_rc" : "");
  const QString preset = safePart(presetCombo_->currentText().trimmed().isEmpty() ? "default" : presetCombo_->currentText().trimmed());
  const QString dirname = qPart + "_vs_" + tPart;

  const QString base = QCoreApplication::applicationDirPath() + "/cache/alignments/" + dirname;
  return base + "/" + dirname + "." + preset + ".paf";
}

void AlignPage::onRunAlign() {
  if (alignRunning_) {
    appendLog("Running: alignment is already in progress.");
    return;
  }
  const QString tf = targetFasta_->text().trimmed();
  const QString qf = queryFasta_->text().trimmed();
  const QString ts = targetSeqCombo_->currentText().trimmed();
  const QString qs = querySeqCombo_->currentText().trimmed();
  if (tf.isEmpty() || qf.isEmpty() || ts.isEmpty() || qs.isEmpty()) {
    QMessageBox::warning(this, "Missing input", "Please select FASTA files and sequences first.");
    return;
  }

  const QString pafPath = computeCachePafPath();
  if (pafPath.isEmpty()) {
    QMessageBox::warning(this, "Path error", "Failed to compute cache PAF path.");
    return;
  }
  cachePathView_->setText(pafPath);

  gapneedle::AlignmentRequest req;
  req.targetFasta = tf.toStdString();
  req.queryFasta = qf.toStdString();
  req.targetSeq = ts.toStdString();
  req.querySeq = qs.toStdString();
  req.outputPafPath = pafPath.toStdString();
  req.preset = presetCombo_->currentText().toStdString();
  req.threads = threads_->value();
  req.reverseTarget = reverseTarget_->isChecked();
  req.reverseQuery = reverseQuery_->isChecked();
  req.reuseExisting = true;

  alignRunning_ = true;
  runBtn_->setEnabled(false);
  appendLog("Running: alignment started in background...");
  emit alignmentStarted(ts, qs);

  auto* watcher = new QFutureWatcher<AlignTaskResult>(this);
  connect(watcher, &QFutureWatcher<AlignTaskResult>::finished, this, [this, watcher, ts, qs, tf, qf]() {
    const AlignTaskResult task = watcher->result();
    watcher->deleteLater();
    alignRunning_ = false;
    runBtn_->setEnabled(true);
    if (!task.configLine.isEmpty()) appendLog(task.configLine);
    if (!task.cacheLine.isEmpty()) appendLog(task.cacheLine);
    if (!task.error.isEmpty()) {
      appendLog(QString("Error: %1").arg(task.error));
      QMessageBox::critical(this, "Alignment failed", task.error);
      emit alignmentFailed(task.error);
      return;
    }
    appendLog(QString("PAF: %1").arg(QString::fromStdString(task.result.pafPath)));
    appendLog(QString("Status: %1").arg(task.result.skipped ? "reused cache" : "newly generated"));
    emit alignmentReady(QString::fromStdString(task.result.pafPath), ts, qs, tf, qf);
  });

  watcher->setFuture(QtConcurrent::run([facade = facade_, req]() {
    AlignTaskResult out;
    try {
      const std::filesystem::path pafPath(req.outputPafPath);
      if (pafPath.has_parent_path()) {
        std::error_code ec;
        std::filesystem::create_directories(pafPath.parent_path(), ec);
        if (ec) {
          out.error = QString("failed to create output directory: %1").arg(QString::fromStdString(ec.message()));
          return out;
        }
      }

      const bool cacheHit = std::filesystem::exists(pafPath);
      out.configLine = QString("Config: target=%1 (%2) | query=%3 (%4) | preset=%5 threads=%6 reverse_target=%7 reverse_query=%8")
                           .arg(QString::fromStdString(req.targetSeq))
                           .arg(QString::fromStdString(req.targetFasta))
                           .arg(QString::fromStdString(req.querySeq))
                           .arg(QString::fromStdString(req.queryFasta))
                           .arg(QString::fromStdString(req.preset))
                           .arg(req.threads)
                           .arg(req.reverseTarget ? "true" : "false")
                           .arg(req.reverseQuery ? "true" : "false");
      if (cacheHit) {
        out.cacheLine = QString("Cache hit: %1 (will reuse if aligner supports reuse)")
                            .arg(QString::fromStdString(req.outputPafPath));
      } else {
        out.cacheLine = QString("Cache miss: %1").arg(QString::fromStdString(req.outputPafPath));
      }

      out.result = facade->align(req);
    } catch (const std::exception& e) {
      out.error = QString::fromUtf8(e.what());
    } catch (...) {
      out.error = "unknown error";
    }
    return out;
  }));
}
