#pragma once

#include "gapneedle/facade.hpp"

#include <QWidget>

#include <vector>

class QLabel;
class QLineEdit;
class QListWidget;
class QPushButton;
class QSpinBox;
class QTextEdit;

class GuidedStitchPage : public QWidget {
  Q_OBJECT

 public:
  explicit GuidedStitchPage(gapneedle::GapNeedleFacade* facade, QWidget* parent = nullptr);

  void setAlignmentContext(const QString& targetFasta,
                           const QString& queryFasta,
                           const QString& targetSeq,
                           const QString& querySeq,
                           const QString& pafPath);

  std::vector<gapneedle::Segment> selectedSegments() const;

 signals:
  void importRequested(bool append);

 private slots:
  void onStartGuide();
  void onChooseCandidate();
  void onBackStep();
  void onStop();
  void onReset();
  void onImport();

 private:
  bool tryApplyClip(gapneedle::GuidedCandidate* candidate, int clipIndex, QString* errorMessage) const;
  void loadSeedCandidates();
  void loadNextCandidates();
  void refreshPathList();
  void refreshCandidateList();
  void setStatus(const QString& text);
  gapneedle::GuidedConstraints constraintsFromUi() const;
  static QString candidateLabel(const gapneedle::GuidedCandidate& c);

 private:
  gapneedle::GapNeedleFacade* facade_{nullptr};

  QLineEdit* pafPath_{nullptr};
  QLabel* contextSummary_{nullptr};
  QLabel* status_{nullptr};

  QSpinBox* nearZeroSpin_{nullptr};
  QSpinBox* maxJumpSpin_{nullptr};
  QSpinBox* minProgressSpin_{nullptr};
  QSpinBox* maxStepsSpin_{nullptr};
  QSpinBox* topKSpin_{nullptr};
  QLineEdit* clipIndexEdit_{nullptr};
  QLabel* clipHint_{nullptr};

  QListWidget* pathList_{nullptr};
  QListWidget* candidateList_{nullptr};
  QTextEdit* detail_{nullptr};
  QPushButton* chooseBtn_{nullptr};

  QString targetFastaPath_;
  QString queryFastaPath_;
  QString targetSeqName_;
  QString querySeqName_;

  std::vector<gapneedle::GuidedCandidate> path_;
  std::vector<gapneedle::GuidedCandidate> currentCandidates_;
  bool stoppedByUser_{false};
};
