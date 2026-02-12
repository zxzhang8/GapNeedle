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
        QString::number(rec.qLen),
        QString::number(rec.qStart),
        QString::number(rec.qEnd),
        QString(QChar(rec.strand)),
        QString::fromStdString(rec.tName),
        QString::number(rec.tLen),
        QString::number(rec.tStart),
        QString::number(rec.tEnd),
        QString::number(rec.matches),
        QString::number(rec.alnLen),
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
  QStringList lines;
  lines << QString("Record: %1 -> %2 (strand %3)").arg(QString::fromStdString(rec.qName), QString::fromStdString(rec.tName), QString(QChar(rec.strand)));
  lines << QString("Query index: %1").arg(r.qPos);
  if (r.qPosOriented.has_value() && r.qPosOriented.value() != r.qPos) {
    lines << QString("Oriented query index: %1").arg(r.qPosOriented.value());
  }
  lines << QString("Query span: %1 - %2").arg(rec.qStart).arg(rec.qEnd);
  lines << QString("Target span: %1 - %2").arg(rec.tStart).arg(rec.tEnd);
  lines << "";
  lines << QString("Reason: %1").arg(QString::fromStdString(r.reason));
  if (r.op != 0) {
    lines << QString("Operation: %1 (len %2, offset %3)").arg(QChar(r.op)).arg(r.opLen).arg(r.opOffset);
  }
  if (r.tPos.has_value()) {
    lines << QString("Mapped target index: %1").arg(r.tPos.value());
  }
  lines << "";
  lines << "Consumed before this index:";
  lines << QString("  Query-consuming: %1").arg(r.qConsumedBefore);
  lines << QString("  Target-consuming: %1").arg(r.tConsumedBefore);
  lines << "";
  lines << "Insertions / deletions before this index:";
  lines << QString("  Insertions (I): %1").arg(countValue(r.countsBefore, 'I'));
  lines << QString("  Deletions (D): %1").arg(countValue(r.countsBefore, 'D'));
  lines << "";
  lines << "CIGAR totals in this record:";
  lines << QString("  Matches (M/= /X): %1").arg(countValue(r.countsTotal, 'M') + countValue(r.countsTotal, '=') + countValue(r.countsTotal, 'X'));
  lines << QString("  Insertions (I): %1").arg(countValue(r.countsTotal, 'I'));
  lines << QString("  Deletions (D): %1").arg(countValue(r.countsTotal, 'D'));
  lines << QString("  Skips (N): %1").arg(countValue(r.countsTotal, 'N'));
  lines << QString("  Soft clips (S): %1").arg(countValue(r.countsTotal, 'S'));
  lines << QString("  Hard clips (H): %1").arg(countValue(r.countsTotal, 'H'));
  lines << QString("  Pads (P): %1").arg(countValue(r.countsTotal, 'P'));
  return lines.join('\n');
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

  mapDetail_->setPlainText(formatMappingDetail(rec, result));
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
