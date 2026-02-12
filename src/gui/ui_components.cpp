#include "ui_components.hpp"

#include "ui_theme.hpp"

#include <QLabel>
#include <QTimer>
#include <QWidget>

#include <algorithm>

namespace gapneedle::ui {

void showToast(QWidget* parent, const QString& message, const QString& level, int durationMs) {
  if (!parent) {
    return;
  }
  auto* toast = new QLabel(message, parent);
  toast->setStyleSheet(QString("background:%1;color:white;padding:8px 12px;border-radius:8px;")
                           .arg(noticeColor(level)));
  toast->setAttribute(Qt::WA_DeleteOnClose, true);
  toast->adjustSize();

  const int x = parent->width() - toast->width() - 24;
  const int y = 18;
  toast->move(std::max(12, x), y);
  toast->show();

  QTimer::singleShot(durationMs, toast, [toast]() { toast->close(); });
}

}  // namespace gapneedle::ui
