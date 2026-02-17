#include "paf_viewer_page.hpp"

#include "gapneedle/mapping_service.hpp"
#include "gapneedle/paf.hpp"

#include <QApplication>
#include <QEvent>
#include <QFormLayout>
#include <QGridLayout>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollBar>
#include <QSpinBox>
#include <QTableWidget>
#include <QTabWidget>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QComboBox>
#include <QClipboard>
#include <QWheelEvent>

#include <algorithm>
#include <random>

namespace {

QString escaped(const std::string& s) {
  return QString::fromStdString(s).toHtmlEscaped();
}

double pct(int value, int total) {
  if (total <= 0) return 0.0;
  const double raw = (100.0 * static_cast<double>(value)) / static_cast<double>(total);
  return std::clamp(raw, 0.0, 100.0);
}

QString statusBadge(const QString& text, const QString& tone) {
  return QString("<span class='badge %1'>%2</span>").arg(tone, text.toHtmlEscaped());
}

QString axisBar(const QString& title,
                int start,
                int end,
                int total,
                const QString& markerLabel,
                std::optional<int> markerPos,
                const QString& fillColor) {
  const double leftPct = pct(start, total);
  const double widthPct = std::max(0.0, pct(end, total) - leftPct);
  const double markerPct = markerPos.has_value() ? pct(markerPos.value(), total) : -1.0;
  QString html;
  html += "<div class='axis-wrap'>";
  html += QString("<div class='axis-title'>%1</div>").arg(title.toHtmlEscaped());
  html += "<div class='track'>";
  html += QString("<div class='range' style='left:%1%%; width:%2%%; background:%3;'></div>")
              .arg(leftPct, 0, 'f', 2)
              .arg(widthPct, 0, 'f', 2)
              .arg(fillColor);
  if (markerPos.has_value()) {
    html += QString("<div class='marker' style='left:%1%%;'></div>").arg(markerPct, 0, 'f', 2);
  }
  html += "</div>";
  html += QString("<div class='axis-meta'>%1-%2 / %3</div>").arg(start).arg(end).arg(total);
  if (markerPos.has_value()) {
    html += QString("<div class='marker-meta'>%1: %2</div>").arg(markerLabel.toHtmlEscaped()).arg(markerPos.value());
  }
  html += "</div>";
  return html;
}

}  // namespace

PafViewerPage::PafViewerPage(QWidget* parent) : QWidget(parent) {
  auto* outer = new QVBoxLayout(this);
  outer->setContentsMargins(12, 12, 12, 12);
  outer->setSpacing(8);

  auto* form = new QFormLayout();
  pafPath_ = new QLineEdit(this);
  auto* loadBtn = new QPushButton("Load PAF", this);
  loadBtn->setObjectName("primaryButton");
  auto* pafRow = new QWidget(this);
  auto* pafRowLayout = new QHBoxLayout(pafRow);
  pafRowLayout->setContentsMargins(0, 0, 0, 0);
  pafRowLayout->setSpacing(6);
  pafRowLayout->addWidget(pafPath_, 1);
  pafRowLayout->addWidget(loadBtn);
  targetSeq_ = new QLineEdit(this);
  querySeq_ = new QLineEdit(this);
  form->addRow("PAF path", pafRow);
  form->addRow("Target sequence", targetSeq_);
  form->addRow("Query sequence", querySeq_);
  outer->addLayout(form);
  connect(loadBtn, &QPushButton::clicked, this, &PafViewerPage::onLoad);

  tabs_ = new QTabWidget(this);

  auto* recordsPage = new QWidget(this);
  auto* recordsLayout = new QVBoxLayout(recordsPage);

  infoLabel_ = new QLabel("No PAF loaded. Run alignment first.", recordsPage);
  infoLabel_->setObjectName("subtitleLabel");
  infoLabel_->setWordWrap(true);
  recordsLayout->addWidget(infoLabel_);

  auto* control = new QGridLayout();
  control->addWidget(new QLabel("Min mapQ", recordsPage), 0, 0);
  mapqSpin_ = new QSpinBox(recordsPage);
  mapqSpin_->setRange(0, 255);
  mapqSpin_->setValue(0);
  control->addWidget(mapqSpin_, 0, 1);

  control->addWidget(new QLabel("Sort by", recordsPage), 1, 0);
  sortCombo_ = new QComboBox(recordsPage);
  sortCombo_->addItems({"qstart", "qend", "tstart", "tend", "matches"});
  control->addWidget(sortCombo_, 1, 1);

  auto* applyBtn = new QPushButton("Apply", recordsPage);
  applyBtn->setObjectName("primaryButton");
  connect(applyBtn, &QPushButton::clicked, this, &PafViewerPage::onApplyFilterSort);
  control->addWidget(applyBtn, 0, 2, 2, 1);

  auto* colorBtn = new QPushButton("Random color", recordsPage);
  connect(colorBtn, &QPushButton::clicked, this, &PafViewerPage::onApplyRandomRowColor);
  control->addWidget(colorBtn, 0, 3, 2, 1);
  control->setColumnStretch(4, 1);
  recordsLayout->addLayout(control);

  table_ = new QTableWidget(recordsPage);
  table_->setObjectName("recordTable");
  table_->setColumnCount(13);
  table_->setHorizontalHeaderLabels(
      {"qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "matches", "aln_len", "mapq", "q-overlap"});
  table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table_->setAlternatingRowColors(true);
  table_->verticalHeader()->setVisible(false);
  table_->horizontalHeader()->setStretchLastSection(true);
  table_->installEventFilter(this);
  table_->viewport()->installEventFilter(this);
  recordsLayout->addWidget(table_, 1);

  auto* mapRow = new QHBoxLayout();
  mapRow->addWidget(new QLabel("Query index", recordsPage));
  qPosSpin_ = new QSpinBox(recordsPage);
  qPosSpin_->setRange(0, 1000000000);
  mapRow->addWidget(qPosSpin_);
  auto* mapBtn = new QPushButton("Map to target", recordsPage);
  mapBtn->setObjectName("primaryButton");
  connect(mapBtn, &QPushButton::clicked, this, &PafViewerPage::onMapQueryPosition);
  mapRow->addWidget(mapBtn);
  mapResultLabel_ = new QLabel("Select a record, then map a query index (requires cg:Z).", recordsPage);
  mapResultLabel_->setWordWrap(true);
  mapRow->addWidget(mapResultLabel_, 1);
  recordsLayout->addLayout(mapRow);

  tabs_->addTab(recordsPage, "Records");

  auto* detailPage = new QWidget(this);
  auto* detailLayout = new QVBoxLayout(detailPage);
  mapDetail_ = new QTextEdit(detailPage);
  mapDetail_->setReadOnly(true);
  mapDetail_->setStyleSheet("QTextEdit{background:#FFFFFF;color:#1D1D1F;border:1px solid #E5E5EA;border-radius:10px;}");
  detailLayout->addWidget(mapDetail_, 1);
  tabs_->addTab(detailPage, "Map details");

  outer->addWidget(tabs_, 1);
}

void PafViewerPage::setContext(const QString& pafPath,
                               const QString& targetSeq,
                               const QString& querySeq,
                               bool autoLoad) {
  pafPath_->setText(pafPath);
  targetSeq_->setText(targetSeq);
  querySeq_->setText(querySeq);
  if (autoLoad) {
    onLoad();
  }
}

void PafViewerPage::onLoad() {
  if (pafPath_->text().trimmed().isEmpty()) {
    QMessageBox::warning(this, "Missing PAF", "Please provide PAF path.");
    return;
  }
  try {
    allRecords_ = gapneedle::parsePaf(pafPath_->text().toStdString(),
                                      targetSeq_->text().toStdString(),
                                      querySeq_->text().toStdString());
    infoLabel_->setText(QString("PAF: %1 · %2 records").arg(pafPath_->text()).arg(allRecords_.size()));
    onApplyFilterSort();
  } catch (const std::exception& e) {
    QMessageBox::critical(this, "Load failed", e.what());
  }
}

void PafViewerPage::onApplyFilterSort() {
  if (allRecords_.empty()) {
    table_->setRowCount(0);
    shownRecords_.clear();
    return;
  }

  shownRecords_.clear();
  const int minMapq = mapqSpin_->value();
  for (const auto& r : allRecords_) {
    if (r.mapq >= minMapq) {
      shownRecords_.push_back(r);
    }
  }

  const QString key = sortCombo_->currentText();
  std::sort(shownRecords_.begin(), shownRecords_.end(), [key](const gapneedle::AlignmentRecord& a, const gapneedle::AlignmentRecord& b) {
    if (key == "qstart") return a.qStart < b.qStart;
    if (key == "qend") return a.qEnd < b.qEnd;
    if (key == "tstart") return a.tStart < b.tStart;
    if (key == "tend") return a.tEnd < b.tEnd;
    if (key == "matches") return a.matches > b.matches;
    return a.qStart < b.qStart;
  });

  populateTable(shownRecords_);
  mapResultLabel_->setText("Select a record, then map a query index (requires cg:Z).");
  mapDetail_->clear();
}

void PafViewerPage::populateTable(const std::vector<gapneedle::AlignmentRecord>& records) {
  const QLocale locale = QLocale::system();
  std::vector<bool> overlaps(records.size(), false);
  for (int i = 0; i < static_cast<int>(records.size()); ++i) {
    for (int j = 0; j < static_cast<int>(records.size()); ++j) {
      if (i == j) {
        continue;
      }
      if (records[i].qStart < records[j].qEnd && records[i].qEnd > records[j].qStart) {
        overlaps[static_cast<std::size_t>(i)] = true;
        break;
      }
    }
  }

  table_->setRowCount(static_cast<int>(records.size()));
  for (int row = 0; row < static_cast<int>(records.size()); ++row) {
    const auto& rec = records[static_cast<std::size_t>(row)];
    QStringList cols = {
        QString::fromStdString(rec.qName),
        locale.toString(rec.qLen),
        locale.toString(rec.qStart),
        locale.toString(rec.qEnd),
        QString(QChar(rec.strand)),
        QString::fromStdString(rec.tName),
        locale.toString(rec.tLen),
        locale.toString(rec.tStart),
        locale.toString(rec.tEnd),
        locale.toString(rec.matches),
        locale.toString(rec.alnLen),
        QString::number(rec.mapq),
        overlaps[static_cast<std::size_t>(row)] ? "Yes" : "No"};

    for (int col = 0; col < cols.size(); ++col) {
      auto* item = new QTableWidgetItem(cols[col]);
      if (col == 1 || col == 2 || col == 3 || col == 6 || col == 7 || col == 8 || col == 9 || col == 10 || col == 11) {
        item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
      }
      table_->setItem(row, col, item);
    }
  }
  table_->resizeColumnsToContents();
}

void PafViewerPage::onApplyRandomRowColor() {
  const int row = table_->currentRow();
  if (row < 0) {
    QMessageBox::information(this, "No selection", "Select a row first.");
    return;
  }

  auto* firstItem = table_->item(row, 0);
  const bool isColored = firstItem && firstItem->data(Qt::UserRole).toBool();
  if (isColored) {
    for (int col = 0; col < table_->columnCount(); ++col) {
      auto* item = table_->item(row, col);
      if (!item) {
        continue;
      }
      item->setBackground(QBrush());
      if (col == 0) {
        item->setData(Qt::UserRole, false);
      }
    }
    return;
  }

  static std::mt19937 rng{std::random_device{}()};
  std::uniform_int_distribution<int> dist(60, 200);
  QColor color(dist(rng), dist(rng), dist(rng), 55);

  for (int col = 0; col < table_->columnCount(); ++col) {
    auto* item = table_->item(row, col);
    if (!item) {
      continue;
    }
    item->setBackground(color);
    if (col == 0) {
      item->setData(Qt::UserRole, true);
    }
  }
}

int PafViewerPage::countValue(const std::unordered_map<char, int>& m, char key) {
  auto it = m.find(key);
  return it == m.end() ? 0 : it->second;
}

QString PafViewerPage::formatMappingDetail(const gapneedle::AlignmentRecord& rec, const gapneedle::MappingResult& r) const {
  const int matchesTotal = countValue(r.countsTotal, 'M') + countValue(r.countsTotal, '=') + countValue(r.countsTotal, 'X');
  const int insertionTotal = countValue(r.countsTotal, 'I');
  const int deletionTotal = countValue(r.countsTotal, 'D');
  const int skipTotal = countValue(r.countsTotal, 'N');
  const int softTotal = countValue(r.countsTotal, 'S');
  const int hardTotal = countValue(r.countsTotal, 'H');
  const int padTotal = countValue(r.countsTotal, 'P');
  const int indelBefore = countValue(r.countsBefore, 'I') + countValue(r.countsBefore, 'D');

  QString reasonText = QString::fromStdString(r.reason);
  QString reasonTone = "warn";
  if (r.reason == "ok") reasonTone = "ok";
  else if (r.reason == "missing_cigar" || r.reason == "bad_cigar") reasonTone = "error";

  QString opText = "N/A";
  if (r.op != 0) {
    opText = QString("%1 (%2 bp, +%3)").arg(QChar(r.op)).arg(r.opLen).arg(r.opOffset);
  }

  QString targetValue = "N/A";
  if (r.tPos.has_value()) {
    targetValue = QString::number(r.tPos.value());
  }

  const QString queryHint = (r.qPosOriented.has_value() && r.qPosOriented.value() != r.qPos)
                                ? QString("oriented: %1").arg(r.qPosOriented.value())
                                : QString("same orientation");

  QString html;
  html += R"(
<html><head><style>
body{font-family:"SF Pro Text","PingFang SC","Segoe UI",sans-serif;background:#FFFFFF !important;color:#1D1D1F !important;margin:0;padding:10px;}
.panel{background:#FFFFFF !important;border:1px solid #E5E5EA;border-radius:12px;padding:14px 16px;}
.title{font-size:16px;font-weight:700;margin:0 0 6px 0;}
.sub{font-size:12px;color:#6E6E73;margin-bottom:10px;}
.badge{display:inline-block;padding:3px 9px;border-radius:999px;font-size:11px;font-weight:700;}
.badge.ok{background:#E9F7EF;color:#1E7A3D;border:1px solid #CBEAD6;}
.badge.warn{background:#FFF5E6;color:#A15A00;border:1px solid #F1D5A8;}
.badge.error{background:#FDEBEC;color:#B3261E;border:1px solid #F4C7CB;}
.section{margin-top:12px;padding:10px;background:#FFFFFF !important;border:1px solid #ECECF0;border-radius:10px;}
.divider{height:1px;background:#4A4A4F;margin:12px 2px 0 2px;opacity:0.85;}
.section-title{font-size:17px;font-weight:800;margin-bottom:9px;letter-spacing:0.2px;}
.section-title.overview{color:#1F5FBF;}
.section-title.coords{color:#0E766E;}
.section-title.spans{color:#7A4A0A;}
.section-title.cigar{color:#6C2E9E;}
.kv{width:100%;border-collapse:collapse;background:#FFFFFF !important;border:1px solid #ECECF0;border-radius:8px;overflow:hidden;}
.kv tr{background:#FFFFFF !important;}
.kv td{font-size:12px;padding:7px 8px;border-bottom:1px solid #F0F0F2;vertical-align:top;background:#FFFFFF !important;color:#1D1D1F !important;}
.kv tr:last-child td{border-bottom:none;}
.kv td:first-child{width:170px;color:#6E6E73 !important;background:#F8F8FA !important;}
.axis-wrap{margin-top:8px;background:#FFFFFF !important;border:1px solid #ECECF0;border-radius:8px;padding:8px;}
.axis-title{font-size:12px;font-weight:650;color:#1D1D1F;margin-bottom:6px;}
.track{position:relative;height:10px;border-radius:999px;background:#ECECF1 !important;overflow:hidden;}
.range{position:absolute;top:0;height:100%;border-radius:999px;}
.marker{position:absolute;top:-3px;width:2px;height:16px;background:#111111;border-radius:2px;}
.axis-meta{font-size:11px;color:#6E6E73;margin-top:6px;}
.marker-meta{font-size:11px;color:#3A3A3C;margin-top:2px;}
</style></head><body>
)";

  html += "<div class='panel'>";
  html += QString("<div class='title'>%1 → %2 <span style='font-size:12px;color:#6E6E73;'>(strand %3)</span></div>")
              .arg(escaped(rec.qName), escaped(rec.tName), QString(QChar(rec.strand)).toHtmlEscaped());
  html += QString("<div class='sub'>%1 &nbsp; reason=%2</div>").arg(statusBadge(reasonText, reasonTone), reasonText.toHtmlEscaped());

  html += "<div class='section'><div class='section-title overview'>Overview</div>";
  html += "<table class='kv'><tbody>";
  html += QString("<tr><td>Query index</td><td><b>%1</b> <span style='color:#6E6E73;'>%2</span></td></tr>")
              .arg(r.qPos)
              .arg(queryHint.toHtmlEscaped());
  html += QString("<tr><td>Mapped target index</td><td><b>%1</b></td></tr>").arg(targetValue.toHtmlEscaped());
  html += QString("<tr><td>Current CIGAR operation</td><td><b>%1</b></td></tr>").arg(opText.toHtmlEscaped());
  html += QString("<tr><td>Consumed before index</td><td>query=%1, target=%2, indel=%3</td></tr>")
              .arg(r.qConsumedBefore)
              .arg(r.tConsumedBefore)
              .arg(indelBefore);
  html += "</tbody></table></div>";

  html += "<div class='divider'></div>";
  html += "<div class='section'><div class='section-title coords'>Coordinates</div>";
  html += axisBar("Query axis", rec.qStart, rec.qEnd, rec.qLen, "selected query", r.qPosOriented, "#64A7FF");
  html += axisBar("Target axis", rec.tStart, rec.tEnd, rec.tLen, "mapped target", r.tPos, "#72D8B2");
  html += "</div>";

  html += "<div class='divider'></div>";
  html += "<div class='section'><div class='section-title spans'>Record Spans</div>";
  html += "<table class='kv'><tbody>";
  html += QString("<tr><td>Query span</td><td>%1 - %2 (len %3)</td></tr>").arg(rec.qStart).arg(rec.qEnd).arg(rec.qLen);
  html += QString("<tr><td>Target span</td><td>%1 - %2 (len %3)</td></tr>").arg(rec.tStart).arg(rec.tEnd).arg(rec.tLen);
  html += QString("<tr><td>Alignment quality</td><td>matches=%1, aln_len=%2, mapq=%3</td></tr>")
              .arg(rec.matches)
              .arg(rec.alnLen)
              .arg(rec.mapq);
  html += "</tbody></table></div>";

  html += "<div class='divider'></div>";
  html += "<div class='section'><div class='section-title cigar'>CIGAR Summary</div>";
  html += "<table class='kv'><tbody>";
  html += QString("<tr><td>Matches (M/= /X)</td><td>%1</td></tr>").arg(matchesTotal);
  html += QString("<tr><td>Insertions (I)</td><td>%1</td></tr>").arg(insertionTotal);
  html += QString("<tr><td>Deletions (D)</td><td>%1</td></tr>").arg(deletionTotal);
  html += QString("<tr><td>Skips (N)</td><td>%1</td></tr>").arg(skipTotal);
  html += QString("<tr><td>Soft clips (S)</td><td>%1</td></tr>").arg(softTotal);
  html += QString("<tr><td>Hard clips (H)</td><td>%1</td></tr>").arg(hardTotal);
  html += QString("<tr><td>Pads (P)</td><td>%1</td></tr>").arg(padTotal);
  html += "</tbody></table></div>";

  html += "</div></body></html>";
  return html;
}

void PafViewerPage::onMapQueryPosition() {
  const int row = table_->currentRow();
  if (row < 0 || row >= static_cast<int>(shownRecords_.size())) {
    QMessageBox::information(this, "No selection", "Select a PAF record first.");
    return;
  }

  const auto& rec = shownRecords_[static_cast<std::size_t>(row)];
  const int qPos = qPosSpin_->value();
  const auto result = gapneedle::mapQueryToTargetDetail(rec, qPos);

  if (result.reason == "missing_cigar") {
    mapResultLabel_->setText("No mapping: PAF record lacks cg:Z (CIGAR).");
    QMessageBox::warning(this, "Missing cg:Z", "This record has no cg:Z, mapping is disabled.");
  } else if (result.reason == "out_of_range") {
    mapResultLabel_->setText("No mapping: query index is outside the record span.");
  } else if (result.reason == "insertion") {
    mapResultLabel_->setText("No mapping: query index falls in an insertion/soft-clip.");
  } else if (!result.tPos.has_value()) {
    mapResultLabel_->setText("No mapping: CIGAR cannot resolve this position.");
  } else {
    mapResultLabel_->setText(
        QString("Mapped: %1[%2] -> %3[%4] (strand %5)")
            .arg(QString::fromStdString(rec.qName))
            .arg(qPos)
            .arg(QString::fromStdString(rec.tName))
            .arg(result.tPos.value())
            .arg(QChar(rec.strand)));
    QApplication::clipboard()->setText(QString::number(result.tPos.value()));
  }

  mapDetail_->setHtml(formatMappingDetail(rec, result));
  tabs_->setCurrentIndex(1);
}

bool PafViewerPage::eventFilter(QObject* watched, QEvent* event) {
  if ((watched == table_ || watched == table_->viewport()) && event->type() == QEvent::Wheel) {
    auto* wheel = static_cast<QWheelEvent*>(event);
    if (wheel->modifiers().testFlag(Qt::ShiftModifier)) {
      auto* h = table_->horizontalScrollBar();
      if (!h) {
        return QWidget::eventFilter(watched, event);
      }

      int delta = 0;
      if (!wheel->pixelDelta().isNull()) {
        delta = wheel->pixelDelta().y();
        if (delta == 0) {
          delta = wheel->pixelDelta().x();
        }
      } else {
        delta = wheel->angleDelta().y();
        if (delta == 0) {
          delta = wheel->angleDelta().x();
        }
      }

      const int step = std::max(1, h->singleStep()) * 3;
      const int direction = (delta < 0) ? 1 : -1;
      h->setValue(h->value() + direction * step);
      return true;
    }
  }
  return QWidget::eventFilter(watched, event);
}
