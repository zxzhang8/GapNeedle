#include "main_window.hpp"

#include "align_page.hpp"
#include "fasta_search_page.hpp"
#include "manual_stitch_page.hpp"
#include "paf_viewer_page.hpp"
#include "ui_components.hpp"

#include <QFrame>
#include <QCloseEvent>
#include <QCoreApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QSplitter>
#include <QStackedWidget>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QWidget>

#include <cstdlib>

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
  setWindowTitle("GapNeedle Qt6");
  resize(1500, 940);
  setupUi();
  setupConnections();
}

void MainWindow::setupUi() {
  auto* central = new QWidget(this);
  auto* root = new QVBoxLayout(central);
  root->setContentsMargins(12, 12, 12, 8);
  root->setSpacing(10);

  auto* header = new QFrame(central);
  header->setObjectName("card");
  auto* headerLayout = new QVBoxLayout(header);
  headerLayout->setContentsMargins(14, 10, 14, 10);
  headerLayout->setSpacing(4);

  headerTitle_ = new QLabel("Align", header);
  headerTitle_->setObjectName("titleLabel");
  headerSubTitle_ = new QLabel("Single-sequence alignment, inspection, and stitch workflow", header);
  headerSubTitle_->setObjectName("subtitleLabel");
  headerLayout->addWidget(headerTitle_);
  headerLayout->addWidget(headerSubTitle_);

  root->addWidget(header);

  auto* splitter = new QSplitter(Qt::Horizontal, central);
  splitter->setChildrenCollapsible(false);

  auto* navCard = new QFrame(splitter);
  navCard->setObjectName("card");
  auto* navLayout = new QVBoxLayout(navCard);
  navLayout->setContentsMargins(8, 10, 8, 10);
  navLayout->setSpacing(8);

  auto* navTitle = new QLabel("Workspace", navCard);
  navTitle->setObjectName("subtitleLabel");
  navLayout->addWidget(navTitle);

  navList_ = new QListWidget(navCard);
  navList_->setObjectName("navList");
  navList_->setFocusPolicy(Qt::NoFocus);
  navList_->addItem("Align");
  navList_->addItem("PAF Viewer");
  navList_->addItem("Manual Stitch");
  navList_->addItem("FASTA Search");
  navList_->setCurrentRow(0);
  navLayout->addWidget(navList_, 1);

  pages_ = new QStackedWidget(splitter);
  alignPage_ = new AlignPage(&facade_, pages_);
  pafViewerPage_ = new PafViewerPage(pages_);
  manualPage_ = new ManualStitchPage(&facade_, pages_);
  searchPage_ = new FastaSearchPage(pages_);

  pages_->addWidget(alignPage_);
  pages_->addWidget(pafViewerPage_);
  pages_->addWidget(manualPage_);
  pages_->addWidget(searchPage_);

  splitter->addWidget(navCard);
  splitter->addWidget(pages_);
  splitter->setStretchFactor(0, 0);
  splitter->setStretchFactor(1, 1);
  splitter->setSizes({240, 1220});

  root->addWidget(splitter, 1);
  setCentralWidget(central);

  statusIcon_ = new QLabel("●", this);
  statusIcon_->setObjectName("statusDot");
  statusIcon_->setFixedWidth(18);
  statusIcon_->setAlignment(Qt::AlignCenter);
  statusBar()->addPermanentWidget(statusIcon_);
  setStatusIcon("idle", "Ready");
}

void MainWindow::setupConnections() {
  connect(navList_, &QListWidget::currentRowChanged, this, &MainWindow::onNavChanged);

  connect(alignPage_,
          &AlignPage::alignmentStarted,
          this,
          [this](const QString& targetSeq, const QString& querySeq) {
            manualPage_->setExternalBusy(true, "alignment is running");
            setStatusIcon("running", QString("Running alignment: %1 <- %2").arg(targetSeq, querySeq));
            statusBar()->showMessage(QString("Running alignment: %1 <- %2").arg(targetSeq, querySeq));
          });

  connect(alignPage_,
          &AlignPage::alignmentReady,
          this,
          [this](const QString& pafPath,
                 const QString& targetSeq,
                 const QString& querySeq,
                 const QString& targetFasta,
                 const QString& queryFasta) {
            pafViewerPage_->setContext(pafPath, targetSeq, querySeq, true);
            manualPage_->setAlignmentContext(targetFasta, queryFasta, targetSeq, querySeq, pafPath);
            manualPage_->setExternalBusy(false);
            setStatusIcon("success", QString("Alignment ready: %1").arg(pafPath));
            statusBar()->showMessage(QString("Alignment ready: %1").arg(pafPath));
            gapneedle::ui::showToast(this, "Alignment completed and context synced", "success");
          });

  connect(alignPage_, &AlignPage::alignmentFailed, this, [this](const QString& errorMessage) {
    manualPage_->setExternalBusy(false);
    setStatusIcon("error", errorMessage);
    statusBar()->showMessage(QString("Alignment failed: %1").arg(errorMessage));
  });

  connect(manualPage_, &ManualStitchPage::checkStarted, this, [this]() {
    alignPage_->setExternalBusy(true, "breakpoint check is running");
    setStatusIcon("running", "Checking breakpoints...");
    statusBar()->showMessage("Checking breakpoints...");
  });

  connect(manualPage_, &ManualStitchPage::checkFinished, this, [this]() {
    alignPage_->setExternalBusy(false);
    setStatusIcon("success", "Breakpoint check completed");
    statusBar()->showMessage("Breakpoint check completed");
  });

  connect(manualPage_, &ManualStitchPage::checkFailed, this, [this](const QString& errorMessage) {
    alignPage_->setExternalBusy(false);
    setStatusIcon("error", errorMessage);
    statusBar()->showMessage(QString("Breakpoint check failed: %1").arg(errorMessage));
  });
}

void MainWindow::onNavChanged(int row) {
  if (row < 0 || row >= pages_->count()) {
    return;
  }
  pages_->setCurrentIndex(row);

  switch (row) {
    case 0:
      headerTitle_->setText("Align");
      headerSubTitle_->setText("Configure and run cached alignment tasks");
      break;
    case 1:
      headerTitle_->setText("PAF Viewer");
      headerSubTitle_->setText("Inspect records, map coordinates, and validate overlap candidates");
      break;
    case 2:
      headerTitle_->setText("Manual Stitch");
      headerSubTitle_->setText("Compose segments, verify breakpoints, and export merged FASTA");
      break;
    case 3:
      headerTitle_->setText("FASTA Search");
      headerSubTitle_->setText("Find query subsequences and inspect hit coordinates");
      break;
    default:
      break;
  }
}

void MainWindow::setStatusIcon(const QString& level, const QString& tooltip) {
  QString color = "#8E8E93";  // idle gray
  if (level == "success") color = "#30D158";
  else if (level == "running") color = "#0A84FF";
  else if (level == "warning") color = "#FF9F0A";
  else if (level == "error") color = "#FF453A";

  statusIcon_->setStyleSheet(QString("color:%1; font-size:14px;").arg(color));
  statusIcon_->setToolTip(tooltip);
}

void MainWindow::closeEvent(QCloseEvent* event) {
  if (alignPage_ && alignPage_->isAlignmentRunning()) {
    statusBar()->showMessage("Alignment is still running. Forcing full shutdown...");
    QCoreApplication::processEvents();
    event->accept();
    std::exit(0);
  }
  QMainWindow::closeEvent(event);
}
