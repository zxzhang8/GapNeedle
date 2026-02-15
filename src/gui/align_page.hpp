#pragma once

#include "gapneedle/facade.hpp"

#include <QHash>
#include <QWidget>

class QCheckBox;
class QComboBox;
class QLineEdit;
class QPushButton;
class QTextEdit;
class QSpinBox;

class AlignPage : public QWidget {
  Q_OBJECT

 public:
  explicit AlignPage(gapneedle::GapNeedleFacade* facade, QWidget* parent = nullptr);
  bool isAlignmentRunning() const { return alignRunning_; }
  void setExternalBusy(bool busy, const QString& reason = QString());

 signals:
  void alignmentStarted(const QString& targetSeq, const QString& querySeq);
  void alignmentReady(const QString& pafPath,
                      const QString& targetSeq,
                      const QString& querySeq,
                      const QString& targetFasta,
                      const QString& queryFasta);
  void alignmentFailed(const QString& errorMessage);

 private slots:
  void onBrowseTarget();
  void onBrowseQuery();
  void onTargetPathEdited();
  void onQueryPathEdited();
  void onRunAlign();

 private:
  void loadSequenceNames(const QString& fastaPath, QComboBox* combo);
  QString computeCachePafPath() const;
  void appendLog(const QString& text);
  QString normalizedPathKey(const QString& path) const;

 private:
  gapneedle::GapNeedleFacade* facade_;
  QLineEdit* targetFasta_{nullptr};
  QLineEdit* queryFasta_{nullptr};
  QComboBox* targetSeqCombo_{nullptr};
  QComboBox* querySeqCombo_{nullptr};
  QComboBox* presetCombo_{nullptr};
  QSpinBox* threads_{nullptr};
  QCheckBox* reverseTarget_{nullptr};
  QCheckBox* reverseQuery_{nullptr};
  QLineEdit* cachePathView_{nullptr};
  QTextEdit* log_{nullptr};
  QPushButton* runBtn_{nullptr};
  bool alignRunning_{false};
  bool externalBusy_{false};
  QString externalBusyReason_;
  QHash<QString, QStringList> fastaNamesCache_;
  int targetLoadToken_{0};
  int queryLoadToken_{0};
};
