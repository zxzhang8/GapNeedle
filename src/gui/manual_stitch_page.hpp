#pragma once

#include "gapneedle/facade.hpp"
#include "gapneedle/fasta_io.hpp"

#include <QMap>
#include <QString>
#include <QStringList>
#include <QWidget>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

class QComboBox;
class QLineEdit;
class QListWidget;
class QPushButton;
class QSpinBox;
class QTextEdit;
class QVBoxLayout;

class ManualStitchPage : public QWidget {
  Q_OBJECT

 public:
  explicit ManualStitchPage(gapneedle::GapNeedleFacade* facade, QWidget* parent = nullptr);

 public slots:
  void setAlignmentContext(const QString& targetFasta,
                           const QString& queryFasta,
                           const QString& targetSeq,
                           const QString& querySeq,
                           const QString& pafPath);

 private slots:
  void onBrowseTarget();
  void onBrowseQuery();
  void onLoadTargetNames();
  void onLoadQueryNames();
  void onAddExtraSource();
  void onSourceChanged();
  void onAddSegment();
  void onRemoveSegment();
  void onMoveSegmentUp();
  void onMoveSegmentDown();
  void onResumeSegment();
  void onSegmentSelectionChanged();
  void onCheckBreakpoints();
  void onExport();
  void onLoadLog();

 private:
  struct SegmentItem {
    QString source;
    QString seqName;
    int start{0};
    int end{0};
    bool reverse{false};
    QString seq;
    QString leftBefore;
    QString leftAfter;
    QString rightBefore;
    QString rightAfter;
  };

  struct ExtraSourceRow {
    QString key;
    QWidget* rowWidget{nullptr};
    QLineEdit* pathEdit{nullptr};
  };

  void appendResult(const QString& text);
  bool loadNamesForSource(const QString& sourceKey, const QString& fastaPath, bool verbose = false);
  QString sourcePath(const QString& sourceKey) const;
  QStringList sourceNames(const QString& sourceKey) const;
  void refreshSourceCombo(const QString& keepKey = QString());
  void refreshSeqCombo();
  void refreshSegments();
  void refreshPreview();
  bool materializeAll(int contextBp);
  bool materializeSegment(SegmentItem& seg, int contextBp);
  QStringList fastaNamesFast(const QString& path) const;
  QString readSegment(const QString& sourceKey, const QString& seqName, int start, int end, bool reverse) const;
  QString junctionPreview(const QString& left, const QString& right, int contextBp) const;
  bool allBreakpointsMatch() const;

 private:
  gapneedle::GapNeedleFacade* facade_;

  QLineEdit* targetFasta_{nullptr};
  QLineEdit* queryFasta_{nullptr};
  QLineEdit* pafPath_{nullptr};
  QSpinBox* contextSpin_{nullptr};

  QVBoxLayout* extraRowsLayout_{nullptr};

  QComboBox* sourceCombo_{nullptr};
  QComboBox* seqCombo_{nullptr};
  QLineEdit* startEdit_{nullptr};
  QLineEdit* endEdit_{nullptr};
  QPushButton* reverseBtn_{nullptr};

  QListWidget* segmentList_{nullptr};
  QTextEdit* preview_{nullptr};
  QTextEdit* detail_{nullptr};
  QTextEdit* result_{nullptr};

  std::vector<SegmentItem> segments_;
  QMap<QString, ExtraSourceRow> extras_;
  QMap<QString, QStringList> namesBySource_;
  mutable std::unordered_map<std::string, gapneedle::FastaMap> fastaCache_;
  int nextExtraId_{1};
};
