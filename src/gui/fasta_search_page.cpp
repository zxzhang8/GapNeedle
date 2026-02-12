#include "fasta_search_page.hpp"

#include "gapneedle/fasta_io.hpp"

#include <QFileDialog>
#include <QFormLayout>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

FastaSearchPage::FastaSearchPage(QWidget* parent) : QWidget(parent) {
  auto* layout = new QVBoxLayout(this);
  layout->setContentsMargins(12, 12, 12, 12);
  layout->setSpacing(8);

  auto* form = new QFormLayout();

  fastaPath_ = new QLineEdit(this);
  auto* fastaRow = new QWidget(this);
  auto* fastaLayout = new QHBoxLayout(fastaRow);
  fastaLayout->setContentsMargins(0, 0, 0, 0);
  auto* browse = new QPushButton("Browse", fastaRow);
  fastaLayout->addWidget(fastaPath_, 1);
  fastaLayout->addWidget(browse);

  query_ = new QLineEdit(this);
  query_->setPlaceholderText("ACGT...");
  form->addRow("FASTA path", fastaRow);
  form->addRow("Query sequence", query_);

  auto* run = new QPushButton("Search", this);
  run->setObjectName("primaryButton");
  auto* summary = new QLabel("No search executed.", this);
  summary->setObjectName("subtitleLabel");

  connect(browse, &QPushButton::clicked, this, [this]() {
    const QString p = QFileDialog::getOpenFileName(this,
                                                    "Select FASTA",
                                                    QString(),
                                                    "FASTA (*.fa *.fasta *.fna);;All files (*)");
    if (!p.isEmpty()) {
      fastaPath_->setText(p);
    }
  });
  connect(run, &QPushButton::clicked, this, &FastaSearchPage::onSearch);

  table_ = new QTableWidget(this);
  table_->setColumnCount(3);
  table_->setHorizontalHeaderLabels({"seq", "start", "end"});
  table_->setAlternatingRowColors(true);
  table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table_->horizontalHeader()->setStretchLastSection(true);

  layout->addLayout(form);
  layout->addWidget(run);
  layout->addWidget(summary);
  layout->addWidget(table_, 1);
}

void FastaSearchPage::onSearch() {
  table_->setRowCount(0);
  const std::string q = query_->text().trimmed().toUpper().toStdString();
  if (q.empty()) {
    QMessageBox::information(this, "Empty query", "Please enter query sequence.");
    return;
  }

  try {
    const auto fasta = gapneedle::readFasta(fastaPath_->text().toStdString());
    int row = 0;
    for (const auto& [name, seq] : fasta) {
      std::size_t pos = 0;
      while (true) {
        pos = seq.find(q, pos);
        if (pos == std::string::npos) {
          break;
        }
        table_->insertRow(row);
        table_->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(name)));
        table_->setItem(row, 1, new QTableWidgetItem(QString::number(static_cast<qulonglong>(pos))));
        table_->setItem(row, 2,
                        new QTableWidgetItem(QString::number(static_cast<qulonglong>(pos + q.size()))));
        ++row;
        ++pos;
      }
    }
  } catch (const std::exception& e) {
    QMessageBox::critical(this, "Search failed", e.what());
  }
}
