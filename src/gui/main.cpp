#include "main_window.hpp"
#include "ui_theme.hpp"

#include <QApplication>
#include <QIcon>

int main(int argc, char* argv[]) {
  QApplication app(argc, argv);
  app.setApplicationName("GapNeedle");
  app.setWindowIcon(QIcon(":/icons/icons/app_logo.svg"));
  gapneedle::ui::applyAppTheme(app);
  MainWindow w;
  w.show();
  return app.exec();
}
