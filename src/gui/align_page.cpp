#include "align_page.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QCompleter>
#include <QCoreApplication>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QTextEdit>
#include <QVBoxLayout>

#include <fstream>

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

}  // namespace

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
  auto* runBtn = new QPushButton("Run alignment", card);
  runBtn->setObjectName("primaryButton");
  auto* clearBtn = new QPushButton("Clear log", card);
  btnRow->addWidget(runBtn);
  btnRow->addWidget(clearBtn);
  btnRow->addStretch(1);
  layout->addLayout(btnRow);

  log_ = new QTextEdit(card);
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
  connect(runBtn, &QPushButton::clicked, this, &AlignPage::onRunAlign);
  connect(clearBtn, &QPushButton::clicked, log_, &QTextEdit::clear);
}

void AlignPage::appendLog(const QString& text) { log_->append(text); }

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
  combo->clear();
  std::ifstream in(fastaPath.toStdString());
  if (!in) {
    return;
  }

  QStringList names;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      QString name = QString::fromStdString(line.substr(1)).trimmed();
      const int sp = name.indexOf(' ');
      if (sp > 0) {
        name = name.left(sp);
      }
      if (!name.isEmpty()) {
        names.push_back(name);
      }
    }
  }
  names.removeDuplicates();
  combo->addItems(names);
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
  QFileInfo pf(pafPath);
  QDir().mkpath(pf.dir().absolutePath());

  const bool cacheHit = pf.exists();
  appendLog(QString("Config: target=%1 (%2) | query=%3 (%4) | preset=%5 threads=%6 reverse_target=%7 reverse_query=%8")
                .arg(ts)
                .arg(tf)
                .arg(qs)
                .arg(qf)
                .arg(presetCombo_->currentText())
                .arg(threads_->value())
                .arg(reverseTarget_->isChecked() ? "true" : "false")
                .arg(reverseQuery_->isChecked() ? "true" : "false"));

  if (cacheHit) {
    appendLog(QString("Cache hit: %1 (will reuse if aligner supports reuse)").arg(pafPath));
  } else {
    appendLog(QString("Cache miss: %1").arg(pafPath));
  }

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

  try {
    const auto result = facade_->align(req);
    appendLog(QString("PAF: %1").arg(QString::fromStdString(result.pafPath)));
    appendLog(QString("Status: %1").arg(result.skipped ? "reused cache" : "newly generated"));
    emit alignmentReady(QString::fromStdString(result.pafPath), ts, qs, tf, qf);
  } catch (const std::exception& e) {
    appendLog(QString("Error: %1").arg(e.what()));
    QMessageBox::critical(this, "Alignment failed", e.what());
  }
}
