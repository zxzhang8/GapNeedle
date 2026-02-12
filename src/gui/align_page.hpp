#pragma once

#include "gapneedle/facade.hpp"

#include <QWidget>

class QCheckBox;
class QComboBox;
class QLineEdit;
class QTextEdit;
class QSpinBox;

class AlignPage : public QWidget {
  Q_OBJECT

 public:
  explicit AlignPage(gapneedle::GapNeedleFacade* facade, QWidget* parent = nullptr);

 signals:
  void alignmentReady(const QString& pafPath,
                      const QString& targetSeq,
                      const QString& querySeq,
                      const QString& targetFasta,
                      const QString& queryFasta);

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
};
