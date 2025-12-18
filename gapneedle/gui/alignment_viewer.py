from __future__ import annotations

from pathlib import Path
from typing import List

from PyQt5 import QtCore, QtWidgets
from qfluentwidgets import InfoBar, SpinBox, PrimaryPushButton, ComboBox, ScrollArea

from gapneedle.stitcher import parse_paf
from gapneedle.gui.theme import WIDGET_FONT_SIZE


class AlignmentViewer(QtWidgets.QWidget):
    """Read-only view of PAF records for the current selection."""

    def __init__(self, parent: QtWidgets.QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("alignmentViewer")

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = ScrollArea(self)
        scroll.setWidgetResizable(True)
        container = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        self.info_label = QtWidgets.QLabel("No PAF loaded. Run alignment first.")
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        control = QtWidgets.QGridLayout()
        control.setHorizontalSpacing(8)
        control.setVerticalSpacing(4)
        control.addWidget(QtWidgets.QLabel("Min mapQ"), 0, 0)
        self.mapq_spin = SpinBox(self)
        self.mapq_spin.setRange(0, 255)
        self.mapq_spin.setValue(0)
        control.addWidget(self.mapq_spin, 0, 1)
        control.addWidget(QtWidgets.QLabel("Sort by"), 1, 0)
        self.sort_combo = ComboBox(self)
        self.sort_combo.addItems(["qstart", "qend", "tstart", "tend", "matches"])
        control.addWidget(self.sort_combo, 1, 1)
        apply_btn = PrimaryPushButton("Apply", self)
        apply_btn.clicked.connect(self._apply_filter_sort)
        control.addWidget(apply_btn, 0, 2, 2, 1)
        control.setColumnStretch(3, 1)
        layout.addLayout(control)

        self.table = QtWidgets.QTableWidget(self)
        self.table.setColumnCount(12)
        self.table.setHorizontalHeaderLabels(
            [
                "qname",
                "qlen",
                "qstart",
                "qend",
                "strand",
                "tname",
                "tlen",
                "tstart",
                "tend",
                "matches",
                "aln_len",
                "mapq",
            ]
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        layout.addWidget(self.table, 1)
        scroll.setWidget(container)
        outer.addWidget(scroll)

        font = self.table.font()
        font.setPointSize(WIDGET_FONT_SIZE)
        self.table.setFont(font)

        self._records: list = []
        self._current_paf: Optional[Path] = None

    def load_paf(self, paf_path: Path, target_seq: str, query_seq: str) -> None:
        if not paf_path or not Path(paf_path).exists():
            InfoBar.info("No PAF", "PAF file not found. Run alignment first.", parent=self)
            return
        records = parse_paf(paf_path, target_seq, query_seq)
        self._records = records
        self._current_paf = paf_path
        self.info_label.setText(f"PAF: {paf_path} · {len(records)} records")
        self._apply_filter_sort()

    def _apply_filter_sort(self) -> None:
        if not self._records:
            self.table.setRowCount(0)
            return
        min_mapq = self.mapq_spin.value()
        filtered = [r for r in self._records if r.mapq >= min_mapq]
        key = self.sort_combo.currentText()
        key_fn = {
            "qstart": lambda r: r.q_start,
            "qend": lambda r: r.q_end,
            "tstart": lambda r: r.t_start,
            "tend": lambda r: r.t_end,
            "matches": lambda r: r.matches,
        }.get(key, lambda r: r.q_start)
        reverse = key == "matches"
        filtered.sort(key=key_fn, reverse=reverse)

        self.table.setRowCount(len(filtered))
        for row, rec in enumerate(filtered):
            values = [
                rec.query,
                f"{rec.q_len:,}",
                f"{rec.q_start:,}",
                f"{rec.q_end:,}",
                rec.strand,
                rec.target,
                f"{rec.t_len:,}",
                f"{rec.t_start:,}",
                f"{rec.t_end:,}",
                f"{rec.matches:,}",
                f"{rec.aln_len:,}",
                rec.mapq,
            ]
            for col, val in enumerate(values):
                item = QtWidgets.QTableWidgetItem(str(val))
                if col in {1, 2, 3, 6, 7, 8, 9, 10, 11}:
                    item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
