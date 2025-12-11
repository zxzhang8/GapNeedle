from __future__ import annotations

from pathlib import Path

from PyQt5 import QtCore, QtGui, QtWidgets

from gapfillet.stitcher import PafRecord
from gapfillet.gui.theme import WIDGET_FONT_SIZE, ICON_SIZE, CONTROL_HEIGHT


class PafIgvViewer(QtWidgets.QWidget):
    """
    Lightweight IGV-like viewer for PAF records.

    Shows target/query tracks, highlights alignment spans, supports zoom/pan,
    and displays basic stats (strand, mapQ, matches/aln_len).
    """

    def __init__(self, parent: QtWidgets.QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("pafIgvViewer")

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        controls = QtWidgets.QHBoxLayout()
        self.info_label = QtWidgets.QLabel("Load a PAF record from the table to view.")
        self.info_label.setWordWrap(True)
        controls.addWidget(self.info_label, 1)
        self.zoom_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.zoom_slider.setRange(-5, 10)  # exp scale
        self.zoom_slider.setValue(0)
        self.zoom_slider.valueChanged.connect(self._update_zoom)
        controls.addWidget(QtWidgets.QLabel("Zoom"))
        controls.addWidget(self.zoom_slider)
        layout.addLayout(controls)

        self.scene = QtWidgets.QGraphicsScene(self)
        self.view = QtWidgets.QGraphicsView(self.scene, self)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing)
        self.view.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        self.view.setMinimumHeight(300)
        layout.addWidget(self.view, 1)

        self.detail = QtWidgets.QTextEdit(self)
        self.detail.setReadOnly(True)
        font = self.detail.font()
        font.setPointSize(WIDGET_FONT_SIZE)
        self.detail.setFont(font)
        self.detail.setMinimumHeight(140)
        layout.addWidget(self.detail)

        self._current_rec: PafRecord | None = None

    def load_record(self, rec: PafRecord) -> None:
        """Render a single PAF record."""
        self._current_rec = rec
        self.info_label.setText(
            f"{rec.query} -> {rec.target} strand={rec.strand} "
            f"mapQ={rec.mapq} matches={rec.matches:,}/{rec.aln_len:,}"
        )
        self._draw(rec)
        self._populate_detail(rec)

    def _populate_detail(self, rec: PafRecord) -> None:
        lines = [
            f"Query: {rec.query} ({rec.q_len:,} bp)",
            f"Target: {rec.target} ({rec.t_len:,} bp)",
            f"Strand: {rec.strand}",
            f"Query span: {rec.q_start:,} - {rec.q_end:,}",
            f"Target span: {rec.t_start:,} - {rec.t_end:,}",
            f"Matches: {rec.matches:,} / aln_len: {rec.aln_len:,} (mapQ {rec.mapq})",
        ]
        self.detail.setPlainText("\n".join(lines))

    def _draw(self, rec: PafRecord) -> None:
        self.scene.clear()
        padding = 40
        track_height = 24
        gap = 40
        max_len = max(rec.t_len, rec.q_len)
        width = max(600, max_len)  # will be scaled down
        scale = (width - 2 * padding) / max_len

        # Target track
        t_y = padding
        self._draw_track("Target", rec.target, rec.t_len, rec.t_start, rec.t_end, t_y, track_height, scale)
        # Query track
        q_y = padding + track_height + gap
        self._draw_track("Query", rec.query, rec.q_len, rec.q_start, rec.q_end, q_y, track_height, scale)

        self.view.fitInView(self.scene.itemsBoundingRect(), QtCore.Qt.KeepAspectRatio)

    def _draw_track(
        self,
        label: str,
        name: str,
        length: int,
        start: int,
        end: int,
        y: int,
        h: int,
        scale: float,
    ) -> None:
        pen = QtGui.QPen(QtGui.QColor("#0f766e"))
        brush = QtGui.QBrush(QtGui.QColor("#7dd3fc"))
        text_pen = QtGui.QPen(QtGui.QColor("#111827"))
        full_rect = QtCore.QRectF(20, y, length * scale, h)
        self.scene.addRect(full_rect, QtGui.QPen(QtGui.QColor("#cbd5e1")))
        aln_rect = QtCore.QRectF(20 + start * scale, y, (end - start) * scale, h)
        self.scene.addRect(aln_rect, pen, brush)
        label_item = self.scene.addText(f"{label}: {name} ({length:,} bp)")
        label_item.setDefaultTextColor(QtGui.QColor("#0f172a"))
        label_item.setPos(full_rect.x(), full_rect.y() - h)

        start_item = self.scene.addText(f"{start:,}")
        start_item.setDefaultTextColor(QtGui.QColor("#475569"))
        start_item.setPos(aln_rect.x(), y + h + 4)
        end_item = self.scene.addText(f"{end:,}")
        end_item.setDefaultTextColor(QtGui.QColor("#475569"))
        end_item.setPos(aln_rect.x() + aln_rect.width(), y + h + 4)

    def _update_zoom(self) -> None:
        value = self.zoom_slider.value()
        factor = 1.2 ** value
        self.view.resetTransform()
        self.view.scale(factor, factor)
        if self._current_rec:
            self.view.centerOn(self.scene.itemsBoundingRect().center())
