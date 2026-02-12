#pragma once

#include <QWidget>

class QLineEdit;
class QTableWidget;

class FastaSearchPage : public QWidget {
  Q_OBJECT

 public:
  explicit FastaSearchPage(QWidget* parent = nullptr);

 private slots:
  void onSearch();

 private:
  QLineEdit* fastaPath_{nullptr};
  QLineEdit* query_{nullptr};
  QTableWidget* table_{nullptr};
};
