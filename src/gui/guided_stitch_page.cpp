#include "guided_stitch_page.hpp"

#include <QFormLayout>
#include <QFrame>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QSpinBox>
#include <QTextEdit>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace {

QString groupTag(const std::string& g) {
  if (g == "strong") return "[STRONG]";
  if (g == "acceptable") return "[OK]";
  if (g == "risk") return "[RISK]";
  return "[CAND]";
}

}  // namespace

GuidedStitchPage::GuidedStitchPage(gapneedle::GapNeedleFacade* facade, QWidget* parent)
    : QWidget(parent), facade_(facade) {
  auto* root = new QVBoxLayout(this);
  root->setContentsMargins(12, 12, 12, 12);
  root->setSpacing(10);

  auto* title = new QLabel("Guided Stitch (PAF-only semi-automation)", this);
  QFont tf = title->font();
  tf.setBold(true);
  title->setFont(tf);
  root->addWidget(title);

  auto* form = new QFormLayout();
  form->setHorizontalSpacing(10);
  form->setVerticalSpacing(8);

  pafPath_ = new QLineEdit(this);
  contextSummary_ = new QLabel("Context: waiting for Align result", this);
  contextSummary_->setWordWrap(true);
  pafPath_->setReadOnly(true);

  form->addRow("Context", contextSummary_);
  form->addRow("PAF path", pafPath_);
  root->addLayout(form);

  auto* cfgRow = new QHBoxLayout();
  nearZeroSpin_ = new QSpinBox(this);
  nearZeroSpin_->setRange(0, 2000000);
  nearZeroSpin_->setValue(1000);
  maxJumpSpin_ = new QSpinBox(this);
  maxJumpSpin_->setRange(0, 500000000);
  maxJumpSpin_->setValue(200000);
  minProgressSpin_ = new QSpinBox(this);
  minProgressSpin_->setRange(1, 10000000);
  minProgressSpin_->setValue(200);
  maxStepsSpin_ = new QSpinBox(this);
  maxStepsSpin_->setRange(1, 2000);
  maxStepsSpin_->setValue(120);
  topKSpin_ = new QSpinBox(this);
  topKSpin_->setRange(1, 100);
  topKSpin_->setValue(12);

  cfgRow->addWidget(new QLabel("Near-0 window", this));
  cfgRow->addWidget(nearZeroSpin_);
  cfgRow->addWidget(new QLabel("Max jump", this));
  cfgRow->addWidget(maxJumpSpin_);
  cfgRow->addWidget(new QLabel("Min progress", this));
  cfgRow->addWidget(minProgressSpin_);
  cfgRow->addWidget(new QLabel("Max steps", this));
  cfgRow->addWidget(maxStepsSpin_);
  cfgRow->addWidget(new QLabel("Top-K", this));
  cfgRow->addWidget(topKSpin_);
  cfgRow->addStretch(1);
  root->addLayout(cfgRow);

  auto* btnRow = new QHBoxLayout();
  auto* startBtn = new QPushButton("Start guide", this);
  startBtn->setObjectName("primaryButton");
  chooseBtn_ = new QPushButton("Choose selected candidate", this);
  chooseBtn_->setObjectName("primaryButton");
  clipIndexEdit_ = new QLineEdit(this);
  clipIndexEdit_->setPlaceholderText("Clip index (optional)");
  clipHint_ = new QLabel("Valid range: select a candidate first", this);
  clipHint_->setObjectName("subtitleLabel");
  auto* backBtn = new QPushButton("Back one step", this);
  auto* stopBtn = new QPushButton("Stop", this);
  auto* resetBtn = new QPushButton("Reset", this);
  auto* importBtn = new QPushButton("Import to Manual Stitch", this);
  importBtn->setObjectName("primaryButton");
  btnRow->addWidget(startBtn);
  btnRow->addWidget(clipIndexEdit_);
  btnRow->addWidget(chooseBtn_);
  btnRow->addWidget(backBtn);
  btnRow->addWidget(stopBtn);
  btnRow->addWidget(resetBtn);
  btnRow->addWidget(importBtn);
  btnRow->addStretch(1);
  root->addLayout(btnRow);
  root->addWidget(clipHint_);

  auto* split = new QHBoxLayout();
  auto* leftCard = new QFrame(this);
  leftCard->setObjectName("card");
  auto* leftLayout = new QVBoxLayout(leftCard);
  auto* leftTitle = new QLabel("Selected path", leftCard);
  leftTitle->setStyleSheet("background:#FFFFFF;");
  leftLayout->addWidget(leftTitle);
  pathList_ = new QListWidget(leftCard);
  leftLayout->addWidget(pathList_, 1);

  auto* midCard = new QFrame(this);
  midCard->setObjectName("card");
  auto* midLayout = new QVBoxLayout(midCard);
  auto* midTitle = new QLabel("Next candidates", midCard);
  midTitle->setStyleSheet("background:#FFFFFF;");
  midLayout->addWidget(midTitle);
  candidateList_ = new QListWidget(midCard);
  midLayout->addWidget(candidateList_, 1);

  auto* rightCard = new QFrame(this);
  rightCard->setObjectName("card");
  auto* rightLayout = new QVBoxLayout(rightCard);
  auto* rightTitle = new QLabel("Details", rightCard);
  rightTitle->setStyleSheet("background:#FFFFFF;");
  rightLayout->addWidget(rightTitle);
  detail_ = new QTextEdit(rightCard);
  detail_->setReadOnly(true);
  rightLayout->addWidget(detail_, 1);

  split->addWidget(leftCard, 1);
  split->addWidget(midCard, 2);
  split->addWidget(rightCard, 2);
  root->addLayout(split, 1);

  status_ = new QLabel("Ready. Run Align first, then start guided stitch.", this);
  status_->setObjectName("subtitleLabel");
  root->addWidget(status_);

  connect(startBtn, &QPushButton::clicked, this, &GuidedStitchPage::onStartGuide);
  connect(chooseBtn_, &QPushButton::clicked, this, &GuidedStitchPage::onChooseCandidate);
  connect(backBtn, &QPushButton::clicked, this, &GuidedStitchPage::onBackStep);
  connect(stopBtn, &QPushButton::clicked, this, &GuidedStitchPage::onStop);
  connect(resetBtn, &QPushButton::clicked, this, &GuidedStitchPage::onReset);
  connect(importBtn, &QPushButton::clicked, this, &GuidedStitchPage::onImport);
  connect(candidateList_, &QListWidget::currentRowChanged, this, [this](int row) {
    if (row < 0 || row >= static_cast<int>(currentCandidates_.size())) {
      detail_->clear();
      clipHint_->setText("Valid range: select a candidate first");
      return;
    }
    const auto& c = currentCandidates_[static_cast<std::size_t>(row)];
    clipHint_->setText(QString("Valid range: (%1, %2]").arg(c.segment.start).arg(c.segment.end));
    std::ostringstream oss;
    oss << "Group: " << c.group << "\n"
        << "Score: " << c.score << "\n"
        << "Segment: " << c.segment.source << ":" << c.segment.seqName << ":" << c.segment.start << ":" << c.segment.end
        << (c.segment.reverse ? ":rc" : "") << "\n"
        << "Target-axis: " << c.axisStart << "-" << c.axisEnd << "\n"
        << "Unclipped axis end: " << c.unclippedAxisEnd << "\n"
        << "Support: " << c.supportCount << "\n"
        << "Record id(s): " << c.recordId << "\n"
        << "Rationale: " << c.rationale << "\n";
    if (c.clippedAt > 0) {
      oss << "Clipped at: " << c.clippedAt << "\n";
    }
    detail_->setPlainText(QString::fromStdString(oss.str()));
  });
}

void GuidedStitchPage::setAlignmentContext(const QString& targetFasta,
                                           const QString& queryFasta,
                                           const QString& targetSeq,
                                           const QString& querySeq,
                                           const QString& pafPath) {
  targetFastaPath_ = targetFasta.trimmed();
  queryFastaPath_ = queryFasta.trimmed();
  targetSeqName_ = targetSeq.trimmed();
  querySeqName_ = querySeq.trimmed();
  pafPath_->setText(pafPath);
  contextSummary_->setText(QString("target=%1 (%2), query=%3 (%4)")
                               .arg(targetSeqName_,
                                    targetFastaPath_,
                                    querySeqName_,
                                    queryFastaPath_));
  setStatus("Alignment context synced. Click 'Start guide' to generate seed candidates.");
}

std::vector<gapneedle::Segment> GuidedStitchPage::selectedSegments() const {
  std::vector<gapneedle::Segment> out;
  out.reserve(path_.size());
  for (const auto& p : path_) {
    out.push_back(p.segment);
  }
  return out;
}

void GuidedStitchPage::onStartGuide() {
  path_.clear();
  currentCandidates_.clear();
  stoppedByUser_ = false;
  refreshPathList();
  loadSeedCandidates();
}

void GuidedStitchPage::onChooseCandidate() {
  const int row = candidateList_->currentRow();
  if (row < 0 || row >= static_cast<int>(currentCandidates_.size())) {
    QMessageBox::information(this, "No selection", "Select one candidate first.");
    return;
  }
  if (static_cast<int>(path_.size()) >= maxStepsSpin_->value()) {
    QMessageBox::warning(this, "Max steps reached", "Configured maximum step count reached.");
    return;
  }
  gapneedle::GuidedCandidate chosen = currentCandidates_[static_cast<std::size_t>(row)];
  const QString clipText = clipIndexEdit_->text().trimmed();
  if (!clipText.isEmpty()) {
    bool ok = false;
    const int clipIndex = clipText.toInt(&ok);
    if (!ok) {
      QMessageBox::warning(this, "Invalid index", "Clip index must be an integer.");
      return;
    }
    QString err;
    if (!tryApplyClip(&chosen, clipIndex, &err)) {
      QMessageBox::warning(this, "Invalid clip index", err);
      return;
    }
  }
  path_.push_back(std::move(chosen));
  clipIndexEdit_->clear();
  refreshPathList();
  loadNextCandidates();
}

void GuidedStitchPage::onBackStep() {
  if (path_.empty()) {
    return;
  }
  path_.pop_back();
  refreshPathList();
  if (path_.empty()) {
    loadSeedCandidates();
  } else {
    loadNextCandidates();
  }
}

void GuidedStitchPage::onStop() {
  stoppedByUser_ = true;
  setStatus("Guide stopped by user. You can import current path or continue.");
}

void GuidedStitchPage::onReset() {
  path_.clear();
  currentCandidates_.clear();
  stoppedByUser_ = false;
  refreshPathList();
  refreshCandidateList();
  detail_->clear();
  clipIndexEdit_->clear();
  clipHint_->setText("Valid range: select a candidate first");
  setStatus("Guide state reset.");
}

void GuidedStitchPage::onImport() {
  if (path_.empty()) {
    QMessageBox::information(this, "Empty path", "No selected segments to import.");
    return;
  }
  emit importRequested(true);
}

void GuidedStitchPage::loadSeedCandidates() {
  const QString targetSeq = targetSeqName_;
  const QString querySeq = querySeqName_;
  const QString pafPath = pafPath_->text().trimmed();
  if (pafPath.isEmpty() || targetSeq.isEmpty() || querySeq.isEmpty()) {
    QMessageBox::warning(this, "Missing context", "Run Align first. Guided Stitch only accepts synced context.");
    return;
  }
  const QFileInfo fi(pafPath);
  if (!fi.exists() || !fi.isFile()) {
    QMessageBox::warning(this, "Missing PAF", QString("PAF file not found:\n%1").arg(pafPath));
    return;
  }
  try {
    gapneedle::GuidedSeedRequest req;
    req.pafPath = pafPath.toStdString();
    req.targetSeq = targetSeq.toStdString();
    req.querySeq = querySeq.toStdString();
    req.maxSeeds = topKSpin_->value();
    req.constraints = constraintsFromUi();
    auto r = facade_->guidedSeed(req);
    currentCandidates_ = std::move(r.candidates);
    refreshCandidateList();
    if (currentCandidates_.empty()) {
      setStatus("No seed candidate found. Try larger near-zero window.");
    } else {
      setStatus(QString("Loaded %1 seed candidate(s). Select one to continue.").arg(currentCandidates_.size()));
    }
  } catch (const std::exception& e) {
    QMessageBox::critical(this, "Guided seed failed", e.what());
  }
}

void GuidedStitchPage::loadNextCandidates() {
  if (path_.empty()) {
    loadSeedCandidates();
    return;
  }
  try {
    const QString pafPath = pafPath_->text().trimmed();
    const QString targetSeq = targetSeqName_;
    const QString querySeq = querySeqName_;
    const QFileInfo fi(pafPath);
    if (pafPath.isEmpty() || !fi.exists() || !fi.isFile()) {
      QMessageBox::warning(this, "Missing PAF", QString("PAF file not found:\n%1").arg(pafPath));
      return;
    }
    gapneedle::GuidedStepRequest req;
    req.pafPath = pafPath.toStdString();
    req.targetSeq = targetSeq.toStdString();
    req.querySeq = querySeq.toStdString();
    req.maxNext = topKSpin_->value();
    req.constraints = constraintsFromUi();
    req.lastAxisEnd = path_.back().axisEnd;
    req.lastChosenIndex = path_.back().clippedAt;
    req.chosenPath.reserve(path_.size());
    for (const auto& c : path_) {
      req.chosenPath.push_back(c.segment);
    }

    auto r = facade_->guidedNext(req);
    currentCandidates_ = std::move(r.candidates);
    refreshCandidateList();
    if (r.exhausted) {
      setStatus("No more monotonic candidate. You reached reachable end; import or backtrack.");
    } else {
      setStatus(QString("Loaded %1 next candidate(s).").arg(currentCandidates_.size()));
    }
  } catch (const std::exception& e) {
    QMessageBox::critical(this, "Guided next failed", e.what());
  }
}

void GuidedStitchPage::refreshPathList() {
  pathList_->clear();
  long long total = 0;
  for (int i = 0; i < static_cast<int>(path_.size()); ++i) {
    const auto& p = path_[i];
    const int addBp = std::max(0, p.axisEnd - p.axisStart);
    total += addBp;
    pathList_->addItem(QString("[%1] %2:%3 %4-%5 | axis %6-%7 | +%8 bp")
                           .arg(i)
                           .arg(QString::fromStdString(p.segment.source))
                           .arg(QString::fromStdString(p.segment.seqName))
                           .arg(p.segment.start)
                           .arg(p.segment.end)
                           .arg(p.axisStart)
                           .arg(p.axisEnd)
                           .arg(addBp));
    if (p.clippedAt > 0) {
      pathList_->item(pathList_->count() - 1)->setToolTip(QString("clipped at %1").arg(p.clippedAt));
    }
  }
  if (!path_.empty()) {
    setStatus(QString("Current path steps: %1, approximate assembled length: %2 bp")
                  .arg(path_.size())
                  .arg(total));
  }
}

void GuidedStitchPage::refreshCandidateList() {
  candidateList_->clear();
  std::sort(currentCandidates_.begin(), currentCandidates_.end(), [](const auto& a, const auto& b) {
    if (a.group != b.group) {
      const auto rank = [](const std::string& g) {
        if (g == "strong") return 0;
        if (g == "acceptable") return 1;
        return 2;
      };
      return rank(a.group) < rank(b.group);
    }
    return a.score > b.score;
  });

  int suffixBaseAxis = -1;
  if (!path_.empty()) {
    suffixBaseAxis = path_.back().axisEnd;
  }

  for (const auto& c : currentCandidates_) {
    QString label = candidateLabel(c);
    auto* item = new QListWidgetItem(label, candidateList_);
    if (suffixBaseAxis >= 0 && c.axisStart == suffixBaseAxis) {
      item->setText(label + QString("  [suffix@%1]").arg(suffixBaseAxis));
      item->setToolTip(QString("Derived suffix candidate from axis %1").arg(suffixBaseAxis));
    }
  }
  chooseBtn_->setEnabled(!currentCandidates_.empty());
  if (currentCandidates_.empty()) {
    clipHint_->setText("Valid range: no candidate");
  }
}

void GuidedStitchPage::setStatus(const QString& text) {
  status_->setText(text);
}

gapneedle::GuidedConstraints GuidedStitchPage::constraintsFromUi() const {
  gapneedle::GuidedConstraints c;
  c.nearZeroWindow = nearZeroSpin_->value();
  c.maxJumpBp = maxJumpSpin_->value();
  c.minProgressBp = minProgressSpin_->value();
  c.maxSteps = maxStepsSpin_->value();
  return c;
}

bool GuidedStitchPage::tryApplyClip(gapneedle::GuidedCandidate* candidate,
                                    int clipIndex,
                                    QString* errorMessage) const {
  if (!candidate) {
    if (errorMessage) *errorMessage = "internal error: null candidate";
    return false;
  }
  const int segStart = candidate->segment.start;
  const int segEnd = candidate->segment.end;
  if (!(clipIndex > segStart && clipIndex <= segEnd)) {
    if (errorMessage) {
      *errorMessage = QString("Index out of range. Expected (%1, %2].").arg(segStart).arg(segEnd);
    }
    return false;
  }

  const int oldAxisEnd = candidate->axisEnd;
  const int oldAxisStart = candidate->axisStart;
  candidate->unclippedAxisEnd = oldAxisEnd;

  int newAxisEnd = oldAxisEnd;
  if (candidate->segment.source == "t") {
    newAxisEnd = clipIndex;
  } else if (candidate->segment.source == "q") {
    const int qSpan = segEnd - segStart;
    const int tSpan = oldAxisEnd - oldAxisStart;
    if (qSpan <= 0 || tSpan <= 0) {
      if (errorMessage) *errorMessage = "Cannot clip this q-candidate due to invalid span.";
      return false;
    }
    const double ratio = static_cast<double>(clipIndex - segStart) / static_cast<double>(qSpan);
    newAxisEnd = oldAxisStart + static_cast<int>(std::lround(ratio * static_cast<double>(tSpan)));
    newAxisEnd = std::max(oldAxisStart + 1, std::min(oldAxisEnd, newAxisEnd));
  } else {
    if (errorMessage) *errorMessage = "Clip is only supported for t/q candidates.";
    return false;
  }

  candidate->segment.end = clipIndex;
  candidate->axisEnd = newAxisEnd;
  candidate->clippedAt = clipIndex;
  candidate->rationale += " | clipped";
  return true;
}

QString GuidedStitchPage::candidateLabel(const gapneedle::GuidedCandidate& c) {
  return QString("%1 score=%2  %3:%4 %5-%6  axis %7-%8  support=%9")
      .arg(groupTag(c.group))
      .arg(QString::number(c.score, 'f', 3))
      .arg(QString::fromStdString(c.segment.source))
      .arg(QString::fromStdString(c.segment.seqName))
      .arg(c.segment.start)
      .arg(c.segment.end)
      .arg(c.axisStart)
      .arg(c.axisEnd)
      .arg(c.supportCount);
}
