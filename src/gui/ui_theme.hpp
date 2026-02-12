#pragma once

#include <QString>

class QApplication;

namespace gapneedle::ui {

struct UiThemeTokens {
  QString fontPrimary{"SF Pro Text"};
  int fontBase{13};
  int fontSmall{12};
  int fontTitle{16};

  QString bgApp{"#F5F5F7"};
  QString bgCard{"#FFFFFF"};
  QString textPrimary{"#1D1D1F"};
  QString textSecondary{"#6E6E73"};
  QString border{"#E5E5EA"};
  QString accent{"#3F678F"};
  QString success{"#30D158"};
  QString warning{"#FF9F0A"};
  QString error{"#FF453A"};

  int radiusCard{10};
  int radiusControl{7};
};

const UiThemeTokens& tokens();
void applyAppTheme(QApplication& app);
QString noticeColor(const QString& level);

}  // namespace gapneedle::ui
