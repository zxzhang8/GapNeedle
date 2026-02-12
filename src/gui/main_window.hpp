#pragma once

#include "gapneedle/facade.hpp"

#include <QMainWindow>

class AlignPage;
class PafViewerPage;
class ManualStitchPage;
class FastaSearchPage;
class QLabel;
class QListWidget;
class QStackedWidget;

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  explicit MainWindow(QWidget* parent = nullptr);

 private slots:
  void onNavChanged(int row);

 private:
  void setupUi();
  void setupConnections();
  void setStatusIcon(const QString& level, const QString& tooltip = QString());

 private:
  gapneedle::GapNeedleFacade facade_;

  QListWidget* navList_{nullptr};
  QStackedWidget* pages_{nullptr};
  QLabel* headerTitle_{nullptr};
  QLabel* headerSubTitle_{nullptr};
  QLabel* statusIcon_{nullptr};

  AlignPage* alignPage_{nullptr};
  PafViewerPage* pafViewerPage_{nullptr};
  ManualStitchPage* manualPage_{nullptr};
  FastaSearchPage* searchPage_{nullptr};
};
