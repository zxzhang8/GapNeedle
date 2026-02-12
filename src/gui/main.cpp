#include "main_window.hpp"
#include "ui_theme.hpp"

#include <QApplication>

int main(int argc, char* argv[]) {
  QApplication app(argc, argv);
  gapneedle::ui::applyAppTheme(app);
  MainWindow w;
  w.show();
  return app.exec();
}
