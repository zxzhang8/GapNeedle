from __future__ import annotations

from pathlib import Path
import random
from typing import List, Optional

from PyQt5 import QtCore, QtGui, QtWidgets
from qfluentwidgets import InfoBar, SpinBox, PrimaryPushButton, ComboBox, ScrollArea

from gapneedle.stitcher import parse_paf, map_query_to_target_detail
from gapneedle.gui.theme import WIDGET_FONT_SIZE


class AlignmentViewer(QtWidgets.QWidget):
    """Read-only view of PAF records for the current selection."""

    def __init__(self, parent: QtWidgets.QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("alignmentViewer")

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self.tabs = QtWidgets.QTabWidget(self)
        outer.addWidget(self.tabs, 1)

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
        color_btn = PrimaryPushButton("Random color", self)
        color_btn.clicked.connect(self._apply_random_row_color)
        control.addWidget(color_btn, 0, 3, 2, 1, alignment=QtCore.Qt.AlignRight)
        control.setColumnStretch(3, 1)
        layout.addLayout(control)

        self.table = QtWidgets.QTableWidget(self)
        self.table.setColumnCount(13)
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
                "q-overlap",
            ]
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        layout.addWidget(self.table, 1)

        mapper = QtWidgets.QHBoxLayout()
        mapper.addWidget(QtWidgets.QLabel("Query index"))
        self.qpos_spin = SpinBox(self)
        self.qpos_spin.setRange(0, 1_000_000_000)
        mapper.addWidget(self.qpos_spin)
        self.map_btn = PrimaryPushButton("Map to target", self)
        self.map_btn.clicked.connect(self._map_query_position)
        mapper.addWidget(self.map_btn)
        self.map_result = QtWidgets.QLabel("Select a record, then map a query index (requires cg:Z).")
        self.map_result.setWordWrap(True)
        mapper.addWidget(self.map_result, 1)
        layout.addLayout(mapper)
        scroll.setWidget(container)
        self.tabs.addTab(scroll, "Records")

        detail_widget = QtWidgets.QWidget()
        detail_layout = QtWidgets.QVBoxLayout(detail_widget)
        detail_layout.setContentsMargins(16, 16, 16, 16)
        detail_layout.setSpacing(10)
        self.map_detail = QtWidgets.QTextEdit(self)
        self.map_detail.setReadOnly(True)
        detail_layout.addWidget(self.map_detail, 1)
        self.tabs.addTab(detail_widget, "Map details")

        font = self.table.font()
        font.setPointSize(WIDGET_FONT_SIZE)
        self.table.setFont(font)
        detail_font = self.map_detail.font()
        detail_font.setPointSize(WIDGET_FONT_SIZE)
        self.map_detail.setFont(detail_font)

        self._records: list = []
        self._shown_records: list = []
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

        self._shown_records = filtered
        overlaps = [False] * len(filtered)
        for i, rec in enumerate(filtered):
            for j, other in enumerate(filtered):
                if i == j:
                    continue
                if rec.q_start < other.q_end and rec.q_end > other.q_start:
                    overlaps[i] = True
                    break
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
                "Yes" if overlaps[row] else "No",
            ]
            for col, val in enumerate(values):
                item = QtWidgets.QTableWidgetItem(str(val))
                if col in {1, 2, 3, 6, 7, 8, 9, 10, 11}:
                    item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
        self.map_result.setText("Select a record, then map a query index (requires cg:Z).")
        self.map_detail.setPlainText("")

    def _apply_random_row_color(self) -> None:
        row = self.table.currentRow()
        if row < 0:
            InfoBar.info("No selection", "Select a row first.", parent=self)
            return
        first_item = self.table.item(row, 0)
        is_colored = bool(first_item and first_item.data(QtCore.Qt.UserRole))
        if is_colored:
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                if item is not None:
                    item.setBackground(QtGui.QBrush())
                    if col == 0:
                        item.setData(QtCore.Qt.UserRole, False)
            return
        color = self._random_row_color()
        for col in range(self.table.columnCount()):
            item = self.table.item(row, col)
            if item is not None:
                item.setBackground(color)
                if col == 0:
                    item.setData(QtCore.Qt.UserRole, True)

    def _map_query_position(self) -> None:
        row = self.table.currentRow()
        if row < 0 or row >= len(self._shown_records):
            InfoBar.info("No selection", "Select a PAF record first.", parent=self)
            return
        rec = self._shown_records[row]
        q_pos = self.qpos_spin.value()
        result = map_query_to_target_detail(rec, q_pos)
        if result.reason == "missing_cigar":
            self.map_result.setText("No mapping: PAF record lacks cg:Z (CIGAR).")
            InfoBar.warning("Missing cg:Z", "This record has no cg:Z, mapping is disabled.", parent=self)
            self.map_detail.setPlainText("Mapping disabled: cg:Z (CIGAR) is missing in this record.")
            return
        if result.reason == "out_of_range":
            self.map_result.setText("No mapping: query index is outside the record span.")
            self.map_detail.setPlainText("No mapping: query index is outside the PAF record span.")
            return
        if result.reason == "insertion":
            self.map_result.setText("No mapping: query index falls in an insertion/soft-clip.")
        elif result.reason != "ok" or result.t_pos is None:
            self.map_result.setText("No mapping: CIGAR cannot resolve this position.")
        else:
            self.map_result.setText(
                f"Mapped: {rec.query}[{q_pos:,}] -> {rec.target}[{result.t_pos:,}] (strand {rec.strand})"
            )
            QtWidgets.QApplication.clipboard().setText(str(result.t_pos))
            InfoBar.success("Copied", "Target index copied to clipboard.", parent=self)
        self.map_detail.setPlainText(self._format_mapping_detail(rec, result))
        self._switch_to_manual_page()

    def showEvent(self, event: QtGui.QShowEvent) -> None:  # type: ignore[override]
        super().showEvent(event)
        self.tabs.setCurrentIndex(0)

    def _switch_to_manual_page(self) -> None:
        main_window = self.window()
        if hasattr(main_window, "manual_page") and hasattr(main_window, "stackedWidget"):
            try:
                main_window.stackedWidget.setCurrentWidget(main_window.manual_page)
            except Exception:
                pass

    def _format_mapping_detail(self, rec, result) -> str:
        total = result.counts_total
        before = result.counts_before
        lines = [
            f"Record: {rec.query} -> {rec.target} (strand {rec.strand})",
            f"Query index: {result.q_pos:,}",
        ]
        if result.q_pos_oriented is not None and result.q_pos_oriented != result.q_pos:
            lines.append(f"Oriented query index: {result.q_pos_oriented:,}")
        lines.append(f"Query span: {rec.q_start:,} - {rec.q_end:,}")
        lines.append(f"Target span: {rec.t_start:,} - {rec.t_end:,}")
        lines.append("")
        lines.append(f"Reason: {result.reason}")
        if result.op:
            lines.append(f"Operation: {result.op} (len {result.op_len}, offset {result.op_offset})")
        if result.t_pos is not None:
            lines.append(f"Mapped target index: {result.t_pos:,}")
        lines.append("")
        lines.append("Consumed before this index:")
        lines.append(f"  Query-consuming: {result.q_consumed_before:,}")
        lines.append(f"  Target-consuming: {result.t_consumed_before:,}")
        lines.append("")
        lines.append("Insertions / deletions before this index:")
        lines.append(f"  Insertions (I): {before.get('I', 0):,}")
        lines.append(f"  Deletions (D): {before.get('D', 0):,}")
        lines.append("")
        lines.append("CIGAR totals in this record:")
        lines.append(f"  Matches (M/= /X): {total.get('M', 0) + total.get('=', 0) + total.get('X', 0):,}")
        lines.append(f"  Insertions (I): {total.get('I', 0):,}")
        lines.append(f"  Deletions (D): {total.get('D', 0):,}")
        lines.append(f"  Skips (N): {total.get('N', 0):,}")
        lines.append(f"  Soft clips (S): {total.get('S', 0):,}")
        lines.append(f"  Hard clips (H): {total.get('H', 0):,}")
        lines.append(f"  Pads (P): {total.get('P', 0):,}")
        return "\n".join(lines)

    def _random_row_color(self) -> QtGui.QColor:
        return QtGui.QColor(
            random.randint(60, 200),
            random.randint(60, 200),
            random.randint(60, 200),
            55,
        )
