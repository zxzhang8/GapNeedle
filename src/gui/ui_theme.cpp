#include "ui_theme.hpp"

#include <QApplication>
#include <QFont>
#include <QFontDatabase>

namespace gapneedle::ui {

const UiThemeTokens& tokens() {
  static UiThemeTokens t;
  return t;
}

QString noticeColor(const QString& level) {
  const auto& t = tokens();
  if (level == "success") return t.success;
  if (level == "warning") return t.warning;
  if (level == "error") return t.error;
  return t.accent;
}

void applyAppTheme(QApplication& app) {
  const auto& t = tokens();
  app.setStyle("Fusion");

  const QStringList families = QFontDatabase::families();
  QString family = "Segoe UI";
  if (families.contains("SF Pro Text")) family = "SF Pro Text";
  else if (families.contains("PingFang SC")) family = "PingFang SC";
  else if (families.contains("Microsoft YaHei UI")) family = "Microsoft YaHei UI";

  QFont base(family, t.fontBase);
  app.setFont(base);

  const QString qss = QString(R"(
QWidget {
  background: %1;
  color: %2;
}
QMainWindow, QDialog {
  background: %1;
}
QFrame#card {
  background: %3;
  border: 1px solid %4;
  border-radius: %5px;
}
QLineEdit, QPlainTextEdit, QTextEdit, QComboBox, QSpinBox, QTableView, QTableWidget, QListWidget, QTabWidget::pane {
  background: %3;
  border: 1px solid %4;
  border-radius: %6px;
}
QLineEdit, QComboBox, QSpinBox {
  min-height: 28px;
  max-height: 28px;
  padding-left: 8px;
}
QTextEdit, QPlainTextEdit {
  padding: 6px 8px;
}
QListWidget {
  padding: 4px 4px;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QTextEdit:focus, QPlainTextEdit:focus {
  border: 1px solid %7;
}
QPushButton {
  background: %3;
  border: 1px solid %4;
  border-radius: %6px;
  min-height: 28px;
  max-height: 28px;
  padding: 0 10px;
}
QPushButton:hover {
  border-color: #C7C7CC;
}
QPushButton#primaryButton {
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #5E7FA6, stop:1 #3F6188);
  color: white;
  border: 1px solid rgba(26, 51, 77, 0.38);
  font-weight: 600;
}
QPushButton#primaryButton:hover {
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #6E8EB3, stop:1 #4A6D95);
  border-color: rgba(30, 55, 83, 0.34);
}
QPushButton#primaryButton:pressed {
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #3A5A80, stop:1 #2E4A6A);
  border-color: rgba(22, 41, 60, 0.42);
}
QPushButton#primaryButton:disabled {
  background: #AAB7C5;
  border-color: rgba(83, 99, 116, 0.28);
  color: #EEF2F6;
}
QLabel#titleLabel {
  font-size: 17px;
  font-weight: 600;
  color: %2;
  background: transparent;
}
QLabel#subtitleLabel {
  font-size: 12px;
  color: %8;
  background: transparent;
}
QComboBox {
  padding-left: 8px;
  padding-right: 28px;
}
QComboBox::drop-down {
  subcontrol-origin: padding;
  subcontrol-position: top right;
  width: 24px;
  border-left: 1px solid %4;
  background: #F2F2F7;
  border-top-right-radius: %6px;
  border-bottom-right-radius: %6px;
}
QComboBox::down-arrow {
  image: url(:/icons/icons/chevron_down.svg);
  width: 10px;
  height: 10px;
}
QComboBox::drop-down:hover {
  background: #EAEAEE;
}
QComboBox::drop-down:pressed {
  background: #E2E2E8;
}
QSpinBox {
  padding-left: 8px;
  padding-right: 28px;
}
QSpinBox::up-button, QSpinBox::down-button {
  subcontrol-origin: border;
  width: 16px;
  border-left: 1px solid %4;
  background: #F2F2F7;
}
QSpinBox::up-button {
  subcontrol-position: top right;
  border-top-right-radius: %6px;
}
QSpinBox::down-button {
  subcontrol-position: bottom right;
  border-top: 1px solid %4;
  border-bottom-right-radius: %6px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
  background: #EAEAEE;
}
QSpinBox::up-button:pressed, QSpinBox::down-button:pressed {
  background: #E2E2E8;
}
QSpinBox::up-arrow {
  image: url(:/icons/icons/chevron_up.svg);
  width: 10px;
  height: 10px;
}
QSpinBox::down-arrow {
  image: url(:/icons/icons/chevron_down.svg);
  width: 10px;
  height: 10px;
}
QStatusBar {
  background: %3;
  border-top: 1px solid %4;
}
QHeaderView::section {
  background: #F2F2F7;
  border: 0;
  border-right: 1px solid %4;
  border-bottom: 1px solid %4;
  padding: 6px;
}
QListWidget#navList {
  background: transparent;
  border: 0;
  outline: 0;
}
QListWidget#navList::item {
  padding: 8px 10px;
  margin: 2px 6px;
  border-radius: 8px;
  border: 0;
  outline: 0;
}
QListWidget#navList::item:selected {
  background: #E7EEF6;
  color: #355D86;
  border: 0;
  outline: 0;
}
QListWidget#navList::item:selected:active,
QListWidget#navList::item:selected:!active,
QListWidget#navList::item:focus {
  border: 0;
  outline: 0;
}
QTableWidget#recordTable QScrollBar:horizontal {
  background: #F2F2F7;
  height: 12px;
  margin: 2px 8px 4px 8px;
  border: 0;
  border-radius: 6px;
}
QTableWidget#recordTable QScrollBar::handle:horizontal {
  background: #C3C7CE;
  min-width: 36px;
  border-radius: 6px;
}
QTableWidget#recordTable QScrollBar::handle:horizontal:hover {
  background: #A8AFB9;
}
QTableWidget#recordTable QScrollBar::add-line:horizontal,
QTableWidget#recordTable QScrollBar::sub-line:horizontal,
QTableWidget#recordTable QScrollBar::add-page:horizontal,
QTableWidget#recordTable QScrollBar::sub-page:horizontal {
  width: 0px;
  background: transparent;
}
)")
                          .arg(t.bgApp, t.textPrimary, t.bgCard, t.border)
                          .arg(t.radiusCard)
                          .arg(t.radiusControl)
                          .arg(t.accent)
                          .arg(t.textSecondary);

  app.setStyleSheet(qss);
}

}  // namespace gapneedle::ui
