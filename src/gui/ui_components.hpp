#pragma once

#include <QString>

class QWidget;

namespace gapneedle::ui {

void showToast(QWidget* parent, const QString& message, const QString& level = "info", int durationMs = 2200);

}  // namespace gapneedle::ui
