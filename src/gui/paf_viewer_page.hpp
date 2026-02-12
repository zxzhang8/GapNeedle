#pragma once

#include "gapneedle/types.hpp"

#include <QWidget>

#include <vector>

class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;
class QTextEdit;
class QTabWidget;
class QComboBox;

class PafViewerPage : public QWidget {
  Q_OBJECT

 public:
  explicit PafViewerPage(QWidget* parent = nullptr);

 public slots:
  void setContext(const QString& pafPath, const QString& targetSeq, const QString& querySeq, bool autoLoad = true);

 private slots:
  void onLoad();
  void onApplyFilterSort();
  void onApplyRandomRowColor();
  void onMapQueryPosition();

 protected:
  bool eventFilter(QObject* watched, QEvent* event) override;

 private:
  void populateTable(const std::vector<gapneedle::AlignmentRecord>& records);
  QString formatMappingDetail(const gapneedle::AlignmentRecord& rec, const gapneedle::MappingResult& r) const;
  static int countValue(const std::unordered_map<char, int>& m, char key);

 private:
  QLineEdit* pafPath_{nullptr};
  QLineEdit* targetSeq_{nullptr};
  QLineEdit* querySeq_{nullptr};

  QLabel* infoLabel_{nullptr};
  QSpinBox* mapqSpin_{nullptr};
  QComboBox* sortCombo_{nullptr};
  QTableWidget* table_{nullptr};

  QSpinBox* qPosSpin_{nullptr};
  QLabel* mapResultLabel_{nullptr};
  QTextEdit* mapDetail_{nullptr};
  QTabWidget* tabs_{nullptr};

  std::vector<gapneedle::AlignmentRecord> allRecords_;
  std::vector<gapneedle::AlignmentRecord> shownRecords_;
};
