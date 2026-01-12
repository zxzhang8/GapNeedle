from __future__ import annotations

import threading
from dataclasses import dataclass, field
import html
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence

from PyQt5 import QtCore, QtGui, QtWidgets
from qfluentwidgets import (
    ComboBox,
    FluentIcon,
    FluentWindow,
    InfoBar,
    LineEdit,
    NavigationItemPosition,
    PrimaryPushButton,
    PushButton,
    ScrollArea,
    SpinBox,
    CheckBox,
    setTheme,
    setThemeColor,
    Theme,
)

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

try:
    import mappy as mp
except ImportError:
    mp = None

from gapneedle.aligner import AlignmentRun, run_minimap2_alignment
from gapneedle.io import read_fasta_sequences, reverse_complement, write_fasta
from gapneedle.gui.theme import (
    ACCENT_FONT_SIZE,
    BASE_FONT_SIZE,
    BUTTON_PADDING,
    CARD_PADDING,
    CONTROL_HEIGHT,
    FONT_CANDIDATES,
    FONT_SCALING,
    HERO_FONT_SIZE,
    ICON_SIZE,
    LIST_PAGE_SIZE,
    DEFAULT_BROWSE_DIR,
    DEFAULT_SAVE_DIR,
    NAV_BUTTON_HEIGHT,
    NAV_EXPAND_WIDTH,
    NAV_ICON_SIZE,
    NAV_FONT_SIZE,
    NAV_MIN_WIDTH,
    PANEL_PADDING,
    UI_SCALING,
    WIDGET_FONT_SIZE,
    WINDOW_HEIGHT,
    WINDOW_WIDTH, INFO_HEIGHT,
)
from gapneedle.gui.alignment_viewer import AlignmentViewer
from qfluentwidgets.components.navigation.navigation_widget import NavigationPushButton


# -----------------------------
# Shared state / helpers
# -----------------------------

@dataclass
class AppState:
    target_fasta: Optional[Path] = None
    query_fasta: Optional[Path] = None
    target_seq: Optional[str] = None
    query_seq: Optional[str] = None
    preset: str = "asm20"
    threads: int = 4
    reverse_query: bool = False
    last_alignment: Optional[AlignmentRun] = None
    name_cache: Dict[Path, List[str]] = field(default_factory=dict)
    index_cache: Dict[Path, Dict[str, "FaiEntry"]] = field(default_factory=dict)

    def remember_alignment(
        self, run: AlignmentRun, preset: str, reverse_query: bool
    ) -> None:
        self.last_alignment = run
        self.target_fasta = Path(run.target_fasta)
        self.query_fasta = Path(run.query_fasta)
        self.target_seq = run.target_seq
        self.query_seq = run.query_seq
        self.preset = preset
        self.reverse_query = reverse_query

    def has_selection(self) -> bool:
        return all([self.target_fasta, self.query_fasta, self.target_seq, self.query_seq])


def _list_fasta_names(path: Path, state: AppState) -> List[str]:
    idx = _load_fai_index(path, state)
    names = list(idx.keys())
    state.name_cache[path] = names
    return names


@dataclass
class FaiEntry:
    length: int
    offset: int
    line_len: int
    line_blen: int
    sequence: Optional[str] = None


def _read_index_from_fai(fai_path: Path) -> Dict[str, FaiEntry]:
    idx: Dict[str, FaiEntry] = {}
    with Path(fai_path).open() as fh:
        for line in fh:
            if not line.strip():
                continue
            name, slen, start, line_len, line_blen, *_ = line.rstrip("\n").split("\t")
            idx[name] = FaiEntry(
                length=int(slen),
                offset=int(start),
                line_len=int(line_len) or int(slen),
                line_blen=int(line_blen) or int(line_len),
            )
    return idx


def _load_fai_index(path: Path, state: AppState) -> Dict[str, FaiEntry]:
    if path in state.index_cache:
        return state.index_cache[path]
    fai_path = Path(str(path) + ".fai")
    if fai_path.exists():
        idx = _read_index_from_fai(fai_path)
    else:
        try:
            _build_fai(path, fai_path)
            idx = _read_index_from_fai(fai_path)
        except Exception:
            if mp is None:
                raise ImportError("mappy is required to read FASTA headers.")
            idx = {
                name.split()[0]: FaiEntry(
                    length=len(seq),
                    offset=0,
                    line_len=len(seq),
                    line_blen=len(seq),
                    sequence=seq,
                )
                for name, seq, _ in mp.fastx_read(str(path), read_comment=False)
            }
    state.index_cache[path] = idx
    state.name_cache[path] = list(idx.keys())
    return idx


def _read_range_from_fasta(
    path: Path, entry: FaiEntry, start: int, end: int, handle: Optional[object] = None
) -> str:
    """
    Efficient subsequence reader using .fai metadata.
    start/end are 0-based, half-open.
    """
    if start >= end:
        return ""
    if entry.sequence is not None:
        return entry.sequence[start:end]
    line_len = entry.line_len or entry.length
    line_blen = entry.line_blen or line_len
    fh = handle or path.open("rb")
    parts: List[bytes] = []
    pos = start
    while pos < end:
        line_off = pos % line_len
        # number of bases we can take until end of current line or requested end
        take = min(line_len - line_off, end - pos)
        offset = entry.offset + (pos // line_len) * line_blen + line_off
        fh.seek(offset)
        chunk = fh.read(take)
        parts.append(chunk)
        pos += take
    if handle is None:
        fh.close()
    return b"".join(parts).decode()


def _build_fai(fasta_path: Path, fai_path: Path) -> None:
    """
    Minimal FASTA index writer (samtools .fai format) to speed up header listing.
    """
    fasta_path = Path(fasta_path)
    fai_path = Path(fai_path)
    with fasta_path.open("rb") as fh, fai_path.open("w") as out:
        seq_name = None
        seq_len = 0
        seq_start = 0
        line_len = None
        line_blen = None
        pos = 0  # running file offset to avoid frequent tell() calls
        for raw_line in fh:
            raw_len = len(raw_line)
            line = raw_line.rstrip(b"\r\n")
            newline_len = raw_len - len(line)
            if line.startswith(b">"):
                if seq_name is not None:
                    out.write(
                        f"{seq_name}\t{seq_len}\t{seq_start}\t{line_len or 0}\t{line_blen or 0}\n"
                    )
                seq_name = line[1:].split(b" ", 1)[0].decode()
                seq_len = 0
                seq_start = pos + raw_len  # start of sequence after header line
                line_len = None
                line_blen = None
            else:
                base_len = len(line)
                seq_len += base_len
                if line_len is None:
                    line_len = base_len
                    line_blen = base_len + (newline_len or 1)
            pos += raw_len
        if seq_name:
            out.write(f"{seq_name}\t{seq_len}\t{seq_start}\t{line_len or 0}\t{line_blen or 0}\n")


def _format_bp(value: int) -> str:
    return f"{value:,} bp"


def _junction_preview_plain(left: str, right: str, context: int, highlight_len: int = 10) -> str:
    ctx = max(0, context)
    l_tail = left[-ctx:]
    r_head = right[:ctx]
    hl = min(highlight_len, len(l_tail), len(r_head))
    l_head_plain = l_tail[: len(l_tail) - hl]
    l_hl = l_tail[len(l_tail) - hl :]
    r_hl = r_head[:hl]
    r_rest = r_head[hl:]
    if hl == 0:
        return f"{l_tail}|{r_head}"
    return f"{l_head_plain}[{l_hl}]|[{r_hl}]{r_rest}"


def _diff_inline(a: str, b: str) -> tuple[str, str]:
    """
    Highlight differences using [] wrappers to avoid ANSI codes in the GUI.
    """
    res_a: List[str] = []
    res_b: List[str] = []
    for x, y in zip(a, b):
        if x == y:
            res_a.append(x)
            res_b.append(y)
        else:
            res_a.append(f"[{x}]")
            res_b.append(f"[{y}]")
    return "".join(res_a), "".join(res_b)


def _diff_html_colored(a: str, b: str) -> tuple[str, str]:
    """
    HTML-highlight differences: matches green, mismatches red; preserves alignment.
    """
    res_a: List[str] = []
    res_b: List[str] = []
    for x, y in zip(a, b):
        xe = html.escape(x)
        ye = html.escape(y)
        if x == y:
            res_a.append(f'<span style="color:green;">{xe}</span>')
            res_b.append(f'<span style="color:green;">{ye}</span>')
        else:
            res_a.append(f'<span style="color:red;">{xe}</span>')
            res_b.append(f'<span style="color:red;">{ye}</span>')
    return "".join(res_a), "".join(res_b)


def _md_wrap(text: str, color: str) -> str:
    return f'<span style="color:{color};font-weight:bold">{text}</span>'


def _html_wrap(text: str, color: str) -> str:
    return f'<span style="color:{color};font-weight:bold">{text}</span>'


def _highlight_diff_md(a: str, b: str) -> tuple[str, str]:
    """Markdown-friendly diff highlight (orange)."""
    res_a: List[str] = []
    res_b: List[str] = []
    for x, y in zip(a, b):
        if x == y:
            res_a.append(x)
            res_b.append(y)
        else:
            res_a.append(_md_wrap(x, "orange"))
            res_b.append(_md_wrap(y, "orange"))
    return "".join(res_a), "".join(res_b)


def _highlight_diff_html(a: str, b: str) -> tuple[str, str]:
    """HTML-friendly diff highlight."""
    res_a: List[str] = []
    res_b: List[str] = []
    for x, y in zip(a, b):
        if x == y:
            res_a.append(x)
            res_b.append(y)
        else:
            res_a.append(_html_wrap(x, "orange"))
            res_b.append(_html_wrap(y, "orange"))
    return "".join(res_a), "".join(res_b)


def _junction_preview_md(left: str, right: str, context: int, highlight_len: int = 10) -> str:
    """
    Markdown preview with colored junction.
    Matching regions are green; mismatches are orange.
    """
    ctx = max(0, context)
    l_tail = left[-ctx:]
    r_head = right[:ctx]
    hl = min(highlight_len, len(l_tail), len(r_head))
    if hl == 0:
        return f"{l_tail}|{r_head}"
    l_head_plain = l_tail[: len(l_tail) - hl]
    l_hl = l_tail[len(l_tail) - hl :]
    r_hl = r_head[:hl]
    r_rest = r_head[hl:]
    color = "green" if l_hl == r_hl else "orange"
    return f"{l_head_plain}{_md_wrap(l_hl, color)}|{_md_wrap(r_hl, color)}{r_rest}"


def _junction_preview_html(left: str, right: str, context: int, highlight_len: int = 10) -> str:
    """
    HTML preview with colored junction.
    """
    ctx = max(0, context)
    l_tail = left[-ctx:]
    r_head = right[:ctx]
    hl = min(highlight_len, len(l_tail), len(r_head))
    if hl == 0:
        return f"{l_tail}|{r_head}"
    l_head_plain = l_tail[: len(l_tail) - hl]
    l_hl = l_tail[len(l_tail) - hl :]
    r_hl = r_head[:hl]
    r_rest = r_head[hl:]
    color = "green" if l_hl == r_hl else "orange"
    return f"{l_head_plain}{_html_wrap(l_hl, color)}|{_html_wrap(r_hl, color)}{r_rest}"


# -----------------------------
# Workers
# -----------------------------


class NameLoader(QtCore.QThread):
    namesReady = QtCore.pyqtSignal(list)
    failed = QtCore.pyqtSignal(str)
    notice = QtCore.pyqtSignal(str)

    def __init__(self, path: Path, state: AppState, parent: Optional[QtCore.QObject] = None):
        super().__init__(parent)
        self.path = path
        self.state = state

    def run(self) -> None:
        try:
            fai_path = Path(str(self.path) + ".fai")
            if not fai_path.exists():
                self.notice.emit("No .fai index found; building one. This may take a moment for large FASTA files.")
            names = _list_fasta_names(self.path, self.state)
            self.namesReady.emit(names)
        except Exception as exc:  # pragma: no cover
            self.failed.emit(str(exc))


class AlignWorker(QtCore.QThread):
    log = QtCore.pyqtSignal(str)
    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(
        self,
        target_fasta: Path,
        query_fasta: Path,
        target_seq: str,
        query_seq: str,
        preset: str,
        threads: int,
        reverse_query: bool,
        parent: Optional[QtCore.QObject] = None,
    ):
        super().__init__(parent)
        self.target_fasta = target_fasta
        self.query_fasta = query_fasta
        self.target_seq = target_seq
        self.query_seq = query_seq
        self.preset = preset
        self.threads = threads
        self.reverse_query = reverse_query

    def run(self) -> None:
        try:
            self.log.emit(
                f"Config: target={self.target_seq} ({self.target_fasta}) | "
                f"query={self.query_seq} ({self.query_fasta}) | "
                f"preset={self.preset} threads={self.threads} reverse={self.reverse_query}"
            )
            run = run_minimap2_alignment(
                target_fasta=self.target_fasta,
                query_fasta=self.query_fasta,
                target_seq=self.target_seq,
                query_seq=self.query_seq,
                preset=self.preset,
                threads=self.threads,
                reverse_query=self.reverse_query,
            )
            if run.skipped:
                self.log.emit(f"Existing PAF detected, skipped run: {run.output_path}")
            else:
                self.log.emit(f"PAF generated: {run.output_path}")
            self.done.emit(run)
        except Exception as exc:  # pragma: no cover
            self.failed.emit(str(exc))


# -----------------------------
# UI widgets
# -----------------------------


class SequencePicker(QtWidgets.QWidget):
    """
    File + sequence picker with search and pagination.
    """

    sequenceChanged = QtCore.pyqtSignal(str)
    fileChanged = QtCore.pyqtSignal(Path)
    namesLoaded = QtCore.pyqtSignal(int)

    def __init__(self, title: str, state: AppState, parent: Optional[QtWidgets.QWidget] = None):
        super().__init__(parent)
        self.state = state
        self.page_size = LIST_PAGE_SIZE
        self.page = 0
        self.names: List[str] = []
        self.filtered: List[str] = []
        self.current_file: Optional[Path] = None
        self.loader: Optional[NameLoader] = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel(title))
        self.file_edit = LineEdit(self)
        row.addWidget(self.file_edit, 1)
        browse = PrimaryPushButton("Browse", self, icon=FluentIcon.FOLDER)
        browse.clicked.connect(self._browse)
        row.addWidget(browse)
        refresh = PushButton("Refresh names", self)
        refresh.clicked.connect(self.load_names)
        row.addWidget(refresh)
        self.refresh_btn = refresh
        self.browse_btn = browse
        layout.addLayout(row)

        search_row = QtWidgets.QHBoxLayout()
        search_row.addWidget(QtWidgets.QLabel("Search / filter"))
        self.search_edit = LineEdit(self)
        self.search_edit.textChanged.connect(self._apply_filter)
        search_row.addWidget(self.search_edit, 1)
        self.page_label = QtWidgets.QLabel("0/0")
        search_row.addWidget(self.page_label)
        prev_btn = PushButton("Prev", self)
        next_btn = PushButton("Next", self)
        prev_btn.clicked.connect(self.prev_page)
        next_btn.clicked.connect(self.next_page)
        search_row.addWidget(prev_btn)
        search_row.addWidget(next_btn)
        layout.addLayout(search_row)

        self.list_view = QtWidgets.QListWidget(self)
        self.list_view.itemSelectionChanged.connect(self._on_pick)
        layout.addWidget(self.list_view, 1)

    def _browse(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select FASTA",
            directory=str(DEFAULT_BROWSE_DIR) if DEFAULT_BROWSE_DIR else "",
            filter="FASTA (*.fa *.fasta *.fna);;All (*)",
        )
        if path:
            self.file_edit.setText(path)
            self.current_file = Path(path)
            self.fileChanged.emit(self.current_file)
            self.load_names()

    def load_names(self) -> None:
        if not self.file_edit.text():
            InfoBar.error("Select file first", "Please choose a FASTA file.", parent=self)
            return
        path = Path(self.file_edit.text())
        self.current_file = path
        self.fileChanged.emit(path)
        # If cached, load instantly
        if path in self.state.name_cache:
            self.names = self.state.name_cache[path]
            self._apply_filter(reset_page=True)
            self.namesLoaded.emit(len(self.names))
            return
        # Async load to avoid blocking UI on large FASTA
        self.refresh_btn.setEnabled(False)
        self.browse_btn.setEnabled(False)
        self.list_view.clear()
        self.page_label.setText("Loading...")
        self.loader = NameLoader(path, self.state)
        self.loader.namesReady.connect(self._on_names_ready)
        self.loader.failed.connect(self._on_names_failed)
        self.loader.notice.connect(self._on_index_notice)
        self.loader.start()

    def _on_names_ready(self, names: List[str]) -> None:
        self.names = names
        self.refresh_btn.setEnabled(True)
        self.browse_btn.setEnabled(True)
        self._apply_filter(reset_page=True)
        self.namesLoaded.emit(len(names))

    def _on_names_failed(self, err: str) -> None:  # pragma: no cover
        self.refresh_btn.setEnabled(True)
        self.browse_btn.setEnabled(True)
        InfoBar.error("Failed to read", err, parent=self)

    def _on_index_notice(self, msg: str) -> None:
        InfoBar.info("Indexing", msg, duration=8000, parent=self)

    def _apply_filter(self, text: str = "", reset_page: bool = False) -> None:
        query = (text or self.search_edit.text()).strip().lower()
        if query:
            self.filtered = [n for n in self.names if query in n.lower()]
        else:
            self.filtered = list(self.names)
        if reset_page:
            self.page = 0
        self._render_page()

    def _render_page(self) -> None:
        total = len(self.filtered)
        page_count = max(1, (total + self.page_size - 1) // self.page_size) if total else 1
        self.page = max(0, min(self.page, page_count - 1))
        start = self.page * self.page_size
        end = start + self.page_size
        slice_names = self.filtered[start:end]
        self.list_view.clear()
        self.list_view.addItems(slice_names)
        self.page_label.setText(f"{self.page + 1}/{page_count} · {total} items")

    def next_page(self) -> None:
        self.page += 1
        self._render_page()

    def prev_page(self) -> None:
        if self.page > 0:
            self.page -= 1
        self._render_page()

    def _on_pick(self) -> None:
        items = self.list_view.selectedItems()
        if not items:
            return
        name = items[0].text()
        self.sequenceChanged.emit(name)


class LogConsole(QtWidgets.QTextEdit):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setMinimumHeight(200)
        font = QtGui.QFont()
        font.setPointSize(WIDGET_FONT_SIZE)
        self.setFont(font)

    def write(self, msg: str) -> None:
        self.append(msg)
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())

    def clear_log(self) -> None:
        self.clear()


class AlignPage(QtWidgets.QWidget):
    def __init__(self, state: AppState, parent: Optional[QtWidgets.QWidget] = None):
        super().__init__(parent)
        self.setObjectName("alignPage")
        self.state = state
        self.worker: Optional[AlignWorker] = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        layout.setSpacing(PANEL_PADDING)

        pickers = QtWidgets.QHBoxLayout()
        self.target_picker = SequencePicker("Target sequence", state)
        self.query_picker = SequencePicker("Query sequence", state)
        self.target_picker.sequenceChanged.connect(self._on_target_seq)
        self.query_picker.sequenceChanged.connect(self._on_query_seq)
        self.target_picker.fileChanged.connect(self._on_target_file)
        self.query_picker.fileChanged.connect(self._on_query_file)
        pickers.addWidget(self.target_picker, 1)
        pickers.addWidget(self.query_picker, 1)
        layout.addLayout(pickers)

        opt_row = QtWidgets.QHBoxLayout()
        self.thread_spin = SpinBox(self)
        self.thread_spin.setRange(1, 64)
        self.thread_spin.setValue(state.threads)
        opt_row.addWidget(QtWidgets.QLabel("Threads"))
        opt_row.addWidget(self.thread_spin)

        self.preset_combo = ComboBox(self)
        self.preset_combo.addItems(["asm5", "asm10", "asm20", "map-ont", "map-pb", "sr", "splice", "ava-ont", "ava-pb"])
        self.preset_combo.setCurrentText(state.preset)
        opt_row.addWidget(QtWidgets.QLabel("Preset"))
        opt_row.addWidget(self.preset_combo)

        self.reverse_check = CheckBox("Reverse-complement query", self)
        self.reverse_check.setChecked(state.reverse_query)
        opt_row.addWidget(self.reverse_check)
        opt_row.addStretch(1)
        layout.addLayout(opt_row)

        btn_row = QtWidgets.QHBoxLayout()
        run_btn = PrimaryPushButton("Run alignment", self, icon=FluentIcon.PLAY)
        run_btn.clicked.connect(self.start_alignment)
        clear_btn = PushButton("Clear log", self, icon=FluentIcon.DELETE)
        clear_btn.clicked.connect(lambda: self.log.clear_log())
        btn_row.addWidget(run_btn)
        btn_row.addWidget(clear_btn)
        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        self.log = LogConsole(self)
        layout.addWidget(self.log, 1)

    def _on_target_seq(self, name: str) -> None:
        self.state.target_seq = name

    def _on_query_seq(self, name: str) -> None:
        self.state.query_seq = name

    def _on_target_file(self, path: Path) -> None:
        self.state.target_fasta = path

    def _on_query_file(self, path: Path) -> None:
        self.state.query_fasta = path

    def start_alignment(self) -> None:
        if self.worker and self.worker.isRunning():
            InfoBar.info("Running", "Alignment already in progress.", parent=self)
            return
        if not self.state.target_fasta or not self.state.query_fasta:
            InfoBar.warning("Files missing", "Select both target and query FASTA.", parent=self)
            return
        if not self.state.target_seq or not self.state.query_seq:
            InfoBar.warning("Sequences missing", "Pick target/query sequences.", parent=self)
            return

        preset = self.preset_combo.currentText() or "asm10"
        threads = self.thread_spin.value()
        reverse_q = self.reverse_check.isChecked()

        self.worker = AlignWorker(
            target_fasta=self.state.target_fasta,
            query_fasta=self.state.query_fasta,
            target_seq=self.state.target_seq,
            query_seq=self.state.query_seq,
            preset=preset,
            threads=threads,
            reverse_query=reverse_q,
        )
        self.worker.log.connect(self.log.write)
        self.worker.failed.connect(lambda err: InfoBar.error("Align failed", err, parent=self))
        self.worker.done.connect(self._align_done)
        self.log.write("Starting alignment...")
        self.worker.start()

    def _align_done(self, run: AlignmentRun) -> None:
        self.state.remember_alignment(run, self.preset_combo.currentText(), self.reverse_check.isChecked())
        status = "Skipped (existing PAF)" if run.skipped else "Done"
        self.log.write(f"{status} · {run.output_path}")
        InfoBar.success("Alignment finished", f"{status}: {run.output_path}", parent=self)
        self._push_to_viewer(run)

    def _push_to_viewer(self, run: AlignmentRun) -> None:
        main_window = self.window()
        viewer = getattr(main_window, "viewer_widget", None)
        if viewer and hasattr(viewer, "load_paf"):
            try:
                viewer.load_paf(run.output_path, run.target_seq, run.query_seq)
            except Exception:
                pass


class ManualStitchPage(QtWidgets.QWidget):
    def __init__(self, state: AppState, parent: Optional[QtWidgets.QWidget] = None):
        super().__init__(parent)
        self.setObjectName("manualStitchPage")
        self.state = state
        self.segments: List[Dict[str, object]] = []
        self._synced_from_state = False
        self.indexes: Dict[str, Dict[str, FaiEntry]] = {"t": {}, "q": {}}
        self.selected_names: Dict[str, Optional[str]] = {"t": None, "q": None}

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        layout.setSpacing(PANEL_PADDING)

        header = QtWidgets.QHBoxLayout()
        header.addWidget(QtWidgets.QLabel("Manual stitch (coordinate mode)"))
        reuse_btn = PrimaryPushButton("Reuse from Align tab", self, icon=FluentIcon.SYNC)
        reuse_btn.clicked.connect(lambda: self._reuse_state(auto_load=False))
        reuse_load_btn = PushButton("Copy & load", self, icon=FluentIcon.DOWNLOAD)
        reuse_load_btn.clicked.connect(lambda: self._reuse_state(auto_load=True))
        header.addWidget(reuse_btn)
        header.addWidget(reuse_load_btn)
        header.addStretch(1)
        layout.addLayout(header)

        form = QtWidgets.QGridLayout()
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(8)

        self.target_path = LineEdit(self)
        self.query_path = LineEdit(self)
        self.context_spin = SpinBox(self)
        self.context_spin.setRange(20, 5000)
        self.context_spin.setValue(200)
        self.target_path.editingFinished.connect(lambda: self._load_index_for_source("t"))
        self.query_path.editingFinished.connect(lambda: self._load_index_for_source("q"))

        form.addWidget(QtWidgets.QLabel("Target FASTA"), 0, 0)
        form.addWidget(self.target_path, 0, 1)
        browse_t = PushButton("Browse", self)
        browse_t.clicked.connect(lambda: self._browse(self.target_path, "t"))
        form.addWidget(browse_t, 0, 2)

        form.addWidget(QtWidgets.QLabel("Query FASTA"), 1, 0)
        form.addWidget(self.query_path, 1, 1)
        browse_q = PushButton("Browse", self)
        browse_q.clicked.connect(lambda: self._browse(self.query_path, "q"))
        form.addWidget(browse_q, 1, 2)

        form.addWidget(QtWidgets.QLabel("Breakpoint context (bp)"), 0, 3)
        form.addWidget(self.context_spin, 0, 4)

        refresh_idx = PrimaryPushButton("Load .fai", self, icon=FluentIcon.DOWNLOAD)
        refresh_idx.clicked.connect(lambda: (self._load_index_for_source("t"), self._load_index_for_source("q")))
        clear_btn = PushButton("Clear segments", self, icon=FluentIcon.DELETE)
        clear_btn.clicked.connect(self._clear_segments)
        form.addWidget(refresh_idx, 2, 3, 1, 2)
        form.addWidget(clear_btn, 3, 3, 1, 2)
        layout.addLayout(form)

        # Segment controls
        segment_row = QtWidgets.QHBoxLayout()
        self.source_combo = ComboBox(self)
        self.source_combo.addItems(["t (target)", "q (query)"])
        self.source_combo.currentIndexChanged.connect(self._on_source_change)
        segment_row.addWidget(self.source_combo)
        self.seq_combo = QtWidgets.QComboBox(self)
        self.seq_combo.setEditable(True)
        self.seq_combo.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.seq_combo.completer().setFilterMode(QtCore.Qt.MatchContains)
        self.seq_combo.completer().setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        self.seq_combo.currentTextChanged.connect(self._on_seq_change)
        segment_row.addWidget(self.seq_combo, 2)
        self.seq_len_label = QtWidgets.QLabel("")
        segment_row.addWidget(self.seq_len_label)
        self.start_edit = LineEdit(self)
        self.start_edit.setPlaceholderText("start")
        self.end_edit = LineEdit(self)
        self.end_edit.setPlaceholderText("end")
        segment_row.addWidget(self.start_edit)
        segment_row.addWidget(self.end_edit)
        self.segment_reverse_check = CheckBox("Reverse-complement", self)
        self.segment_reverse_check.setToolTip("Apply reverse-complement to this segment (query only).")
        self.segment_reverse_check.setEnabled(False)
        segment_row.addWidget(self.segment_reverse_check)
        add_btn = PrimaryPushButton("Add segment", self, icon=FluentIcon.ADD)
        add_btn.clicked.connect(self._add_segment)
        segment_row.addWidget(add_btn)
        segment_row.addStretch(1)
        layout.addLayout(segment_row)

        self.segment_list = QtWidgets.QListWidget(self)
        self.segment_list.itemSelectionChanged.connect(self._on_segment_select)
        layout.addWidget(self.segment_list, 1)

        detail_split = QtWidgets.QHBoxLayout()
        self.preview = QtWidgets.QTextEdit(self)
        self.preview.setReadOnly(True)
        self.preview.setMinimumHeight(INFO_HEIGHT)
        font = QtGui.QFont()
        font.setPointSize(WIDGET_FONT_SIZE)
        self.preview.setFont(font)
        detail_split.addWidget(self.preview, 2)

        self.segment_detail = QtWidgets.QTextEdit(self)
        self.segment_detail.setReadOnly(True)
        self.segment_detail.setMinimumWidth(320)
        self.segment_detail.setFont(font)
        self.segment_detail.setPlaceholderText("Select a segment to view its context.")
        detail_split.addWidget(self.segment_detail, 1)

        layout.addLayout(detail_split)

        footer = QtWidgets.QHBoxLayout()
        check_btn = PushButton("Check breakpoints", self, icon=FluentIcon.CHECKBOX)
        check_btn.clicked.connect(self._check_breakpoints)
        footer.addWidget(check_btn)
        export_btn = PrimaryPushButton("Export merged FASTA", self, icon=FluentIcon.SAVE)
        export_btn.clicked.connect(self._export)
        footer.addWidget(export_btn)
        move_up = PushButton("Move up", self, icon=FluentIcon.LEFT_ARROW)
        move_up.clicked.connect(lambda: self._move_segment(-1))
        move_down = PushButton("Move down", self, icon=FluentIcon.ARROW_DOWN)
        move_down.clicked.connect(lambda: self._move_segment(1))
        footer.addWidget(move_up)
        footer.addWidget(move_down)
        remove_btn = PushButton("Remove selected", self, icon=FluentIcon.DELETE)
        remove_btn.clicked.connect(self._remove_selected)
        footer.addWidget(remove_btn)
        self.length_label = QtWidgets.QLabel("")
        footer.addWidget(self.length_label)
        footer.addStretch(1)
        layout.addLayout(footer)

    def showEvent(self, event: QtGui.QShowEvent) -> None:  # type: ignore[override]
        super().showEvent(event)
        # Auto-fill once from align state if empty.
        if not self._synced_from_state and self.state.has_selection():
            self._reuse_state(auto_load=False, silent=True)
            self._synced_from_state = True

    def _browse(self, line_edit: LineEdit, source: str) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select FASTA",
            directory=str(DEFAULT_BROWSE_DIR) if DEFAULT_BROWSE_DIR else "",
            filter="FASTA (*.fa *.fasta *.fna);;All (*)",
        )
        if path:
            line_edit.setText(path)
            self._load_index_for_source(source)

    def _reuse_state(self, auto_load: bool = False, silent: bool = False) -> None:
        if not self.state.has_selection():
            InfoBar.info("No selection", "Run Align tab first.", parent=self)
            return
        self.target_path.setText(str(self.state.target_fasta))
        self.query_path.setText(str(self.state.query_fasta))
        if self.state.target_seq:
            self.selected_names["t"] = self.state.target_seq
        if self.state.query_seq:
            self.selected_names["q"] = self.state.query_seq
        # Preload indexes for quick picking.
        if auto_load:
            self._load_index_for_source("t", info=not silent)
            self._load_index_for_source("q", info=not silent)
            self._populate_seq_options()
            for key, name in self.selected_names.items():
                if name:
                    self.seq_combo.setCurrentText(name)
        elif not silent:
            InfoBar.success("Filled", "Copied from Align tab.", parent=self)

    def _get_path_for_source(self, src: str) -> Optional[Path]:
        text = self.target_path.text() if src == "t" else self.query_path.text()
        if not text:
            return None
        return Path(text)

    def _load_index_for_source(self, src: str, info: bool = True) -> None:
        path = self._get_path_for_source(src)
        if path is None or not path.exists():
            if info:
                InfoBar.error("Missing file", f"Set {'target' if src == 't' else 'query'} FASTA first.", parent=self)
            return
        if info and not Path(str(path) + ".fai").exists():
            InfoBar.info(
                "Indexing",
                "No .fai index found; building one. This may take a moment.",
                duration=8000,
                parent=self,
            )
        try:
            idx = _load_fai_index(path, self.state)
            self.indexes[src] = idx
            if self._current_source_key() == src:
                self._populate_seq_options()
            if info:
                InfoBar.success(
                    "Index ready",
                    f"{'Target' if src == 't' else 'Query'}: {len(idx):,} sequences loaded.",
                    parent=self,
                )
        except Exception as exc:  # pragma: no cover
            if info:
                InfoBar.error("Index failed", str(exc), parent=self)

    def _current_source_key(self) -> str:
        return "t" if self.source_combo.currentIndex() == 0 else "q"

    def _on_source_change(self) -> None:
        src = self._current_source_key()
        self.segment_reverse_check.setEnabled(src == "q")
        self._populate_seq_options()

    def _on_seq_change(self, name: str) -> None:
        src = self._current_source_key()
        self.selected_names[src] = name or None
        idx = self.indexes.get(src, {})
        if name and name in idx:
            self.seq_len_label.setText(_format_bp(idx[name].length))
        else:
            self.seq_len_label.setText("")

    def _populate_seq_options(self) -> None:
        src = self._current_source_key()
        idx = self.indexes.get(src, {})
        names = list(idx.keys())
        current = self.seq_combo.currentText() or (self.selected_names.get(src) or "")
        self.seq_combo.blockSignals(True)
        self.seq_combo.clear()
        self.seq_combo.addItems(names)
        if current in names:
            self.seq_combo.setCurrentText(current)
        elif names:
            self.seq_combo.setCurrentIndex(0)
        self.seq_combo.blockSignals(False)
        # update completer model to allow contains search
        self.seq_combo.completer().setModel(self.seq_combo.model())
        if current:
            self.seq_combo.setCurrentText(current)
        self._on_seq_change(self.seq_combo.currentText())

    def _add_segment(self) -> None:
        src_key = self._current_source_key()
        path = self._get_path_for_source(src_key)
        if path is None or not path.exists():
            InfoBar.warning("Missing file", "Select FASTA and load index first.", parent=self)
            return
        if not self.indexes.get(src_key):
            self._load_index_for_source(src_key)
            if not self.indexes.get(src_key):
                return
        try:
            start = int(self.start_edit.text())
            end = int(self.end_edit.text())
        except ValueError:
            InfoBar.warning("Invalid", "start/end must be integers.", parent=self)
            return
        seq_name = self.seq_combo.currentText().strip()
        idx = self.indexes.get(src_key, {})
        if not seq_name or seq_name not in idx:
            InfoBar.warning("Sequence", "Pick a valid sequence name.", parent=self)
            return
        seq_len = idx[seq_name].length
        if not (0 <= start < end <= seq_len):
            InfoBar.warning("Out of range", f"Interval invalid for seq length {seq_len:,}", parent=self)
            return
        ctx = self.context_spin.value()
        segment = {
            "src": src_key,
            "name": seq_name,
            "start": start,
            "end": end,
            "length": seq_len,
            "context": ctx,
            "reverse": bool(self.segment_reverse_check.isChecked()) if src_key == "q" else False,
        }
        self.selected_names[src_key] = seq_name
        insert_at = self.segment_list.currentRow()
        if insert_at < 0 or insert_at >= len(self.segments):
            self.segments.append(segment)
            insert_at = len(self.segments) - 1
        else:
            self.segments.insert(insert_at + 1, segment)
            insert_at = insert_at + 1
        self.start_edit.clear()
        self.end_edit.clear()
        self._refresh_segments()
        self.segment_list.setCurrentRow(insert_at)

    def _materialize_segment(
        self, seg: Dict[str, object], ctx: int, handles: Dict[Path, object]
    ) -> Dict[str, object]:
        src = seg["src"]  # type: ignore[index]
        assert src in {"t", "q"}
        path = self._get_path_for_source(src)  # type: ignore[arg-type]
        if path is None:
            raise FileNotFoundError("FASTA path missing.")
        idx = self.indexes.get(src, {})
        name = seg.get("name")
        if name not in idx:
            raise KeyError(f"Sequence {name} not in index.")
        entry = idx[name]  # type: ignore[index]
        fh = handles.get(path)
        if entry.sequence is None:
            if fh is None:
                fh = path.open("rb")
                handles[path] = fh
        else:
            fh = None

        start = int(seg["start"])
        end = int(seg["end"])
        seq_len = entry.length
        reverse = bool(seg.get("reverse")) and src == "q"
        if reverse:
            # Interpret start/end on the reverse-complemented sequence.
            def rc_slice(a: int, b: int) -> str:
                a = max(0, a)
                b = min(seq_len, b)
                if a >= b:
                    return ""
                orig_start = seq_len - b
                orig_end = seq_len - a
                return reverse_complement(
                    _read_range_from_fasta(path, entry, orig_start, orig_end, handle=fh)
                )

            seq = rc_slice(start, end)
            left_before = rc_slice(max(0, start - ctx), start)
            left_after = rc_slice(start, min(seq_len, start + ctx))
            right_before = rc_slice(max(0, end - ctx), end)
            right_after = rc_slice(end, min(seq_len, end + ctx))
        else:
            seq = _read_range_from_fasta(path, entry, start, end, handle=fh)
            left_before = _read_range_from_fasta(path, entry, max(0, start - ctx), start, handle=fh)
            left_after = _read_range_from_fasta(
                path, entry, start, min(seq_len, start + ctx), handle=fh
            )
            right_before = _read_range_from_fasta(path, entry, max(0, end - ctx), end, handle=fh)
            right_after = _read_range_from_fasta(path, entry, end, min(seq_len, end + ctx), handle=fh)

        enriched = dict(seg)
        enriched.update(
            {
                "seq": seq,
                "context": ctx,
                "left_before": left_before,
                "left_after": left_after,
                "right_before": right_before,
                "right_after": right_after,
            }
        )
        return enriched

    def _first_name_from_segments(self, src: str) -> Optional[str]:
        for seg in self.segments:
            if seg.get("src") == src and seg.get("name"):
                return str(seg["name"])
        return None

    def _format_segment_label(self, seg: Dict[str, object]) -> str:
        name = seg.get("name", "?")
        rc_flag = " (RC)" if seg.get("reverse") else ""
        return f"{seg.get('src')}:{name}{rc_flag} {seg.get('start')}-{seg.get('end')}"

    def _refresh_segments(self) -> None:
        self.segment_list.clear()
        total = 0
        for idx, seg in enumerate(self.segments):
            span_len = seg["end"] - seg["start"]
            name = seg.get("name", "?")
            rc_flag = " (RC)" if seg.get("reverse") else ""
            self.segment_list.addItem(
                f"[{idx}] {seg['src']}:{name}{rc_flag} {seg['start']}-{seg['end']} · {span_len:,}bp"
            )
            total += span_len
        self.length_label.setText(f"{len(self.segments)} segments · {_format_bp(total)}")
        self._update_preview()
        if self.segments:
            self.segment_list.setCurrentRow(len(self.segments) - 1)
        else:
            self.segment_detail.clear()

    def _update_preview(self) -> None:
        self.preview.clear()
        if len(self.segments) < 2:
            self.preview.setText("Add two or more segments to show breakpoint context.")
            return
        if not all("left_before" in seg for seg in self.segments):
            lines = ["Sequence data will be read at export time. Current segments:"]
            for seg in self.segments:
                lines.append(
                    f"- {self._format_segment_label(seg)}"
                )
            self.preview.setText("\n".join(lines))
            return
        for i in range(len(self.segments) - 1):
            left = self.segments[i]
            right = self.segments[i + 1]
            ctx_val = left.get("context", self.context_spin.value())  # type: ignore[arg-type]
            preview = _junction_preview_plain(left["right_before"], right["left_after"], ctx_val)
            self.preview.append(
                f"[{i}] {self._format_segment_label(left)} -> {self._format_segment_label(right)}\n{preview}"
            )
            left_slice = left["right_before"][-min(50, len(left["right_before"])) :]
            right_slice = right["left_before"][-min(50, len(right["left_before"])) :]
            if left_slice and right_slice:
                if left_slice == right_slice:
                    self.preview.append('<span style="color:green;">Left flanks match ✓</span>')
                else:
                    hi_prev, hi_curr = _diff_html_colored(left_slice, right_slice)
                    self.preview.append(
                        f'<span style="color:red;">Left flanks differ:</span><br>'
                        f'{hi_prev}<br>{hi_curr}'
                    )
            r_a = left["right_after"][: min(50, len(left["right_after"]))]
            r_b = right["left_after"][: min(50, len(right["left_after"]))]
            if r_a and r_b:
                if r_a == r_b:
                    self.preview.append('<span style="color:green;">Right flanks match ✓</span>\n')
                else:
                    hi_prev, hi_curr = _diff_html_colored(r_a, r_b)
                    self.preview.append(
                        f'<span style="color:red;">Right flanks differ:</span><br>'
                        f'{hi_prev}<br>{hi_curr}<br>'
                    )
        self._on_segment_select()

    def _clear_segments(self) -> None:
        self.segments.clear()
        self._refresh_segments()

    def _remove_selected(self) -> None:
        row = self.segment_list.currentRow()
        if row < 0 or row >= len(self.segments):
            return
        del self.segments[row]
        self._refresh_segments()

    def _move_segment(self, delta: int) -> None:
        row = self.segment_list.currentRow()
        if row < 0 or row >= len(self.segments):
            return
        new_row = row + delta
        if not (0 <= new_row < len(self.segments)):
            return
        self.segments[row], self.segments[new_row] = self.segments[new_row], self.segments[row]
        self._refresh_segments()
        self.segment_list.setCurrentRow(new_row)

    def _on_segment_select(self) -> None:
        row = self.segment_list.currentRow()
        if row < 0 or row >= len(self.segments):
            self.segment_detail.clear()
            return
        seg = self.segments[row]
        ctx = seg.get("context", self.context_spin.value())  # type: ignore[arg-type]
        span_len = seg["end"] - seg["start"]
        lines: List[str] = [
            f"Segment [{row}] {self._format_segment_label(seg)} ({span_len:,} bp)"
        ]
        if "seq" in seg:
            lines += [
                "",
                f"Left context ({ctx} bp before/after start):",
                seg.get("left_before", ""),
                seg.get("left_after", ""),
                "",
                f"Right context ({ctx} bp before/after end):",
                seg.get("right_before", ""),
                seg.get("right_after", ""),
            ]
        else:
            lines.append("")
            lines.append("Sequence not loaded yet; content will be fetched on export.")
        self.segment_detail.setPlainText("\n".join(lines))

    def _export(self) -> None:
        if not self.segments:
            InfoBar.info("No segments", "Add at least one segment.", parent=self)
            return
        ctx = self.context_spin.value()
        enriched = self._materialize_all(ctx)
        if enriched is None:
            return

        merged = "".join(seg["seq"] for seg in enriched)
        target_name = self.selected_names.get("t") or self._first_name_from_segments("t") or "target"
        query_name = self.selected_names.get("q") or self._first_name_from_segments("q") or "query"
        suggested_name = f"{target_name}+{query_name}"
        base_dir = DEFAULT_SAVE_DIR or DEFAULT_BROWSE_DIR
        suggested_path = Path(base_dir) / f"{suggested_name}.fa" if base_dir else f"{suggested_name}.fa"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save merged FASTA",
            str(suggested_path),
            "FASTA (*.fa *.fasta);;All (*)",
        )
        if not path:
            return
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta({suggested_name: merged}, out_path)

        log_path = out_path.with_suffix(".md")
        log_lines = self._build_log_lines(
            enriched=enriched,
            merged_len=len(merged),
            out_path=out_path,
            suggested_name=suggested_name,
            target_name=target_name,
            query_name=query_name,
            ctx=ctx,
            as_html=False,
        )
        log_path.write_text("\n".join(log_lines))

        # Persist enriched sequences for preview after export without re-reading.
        self.segments = enriched
        self._refresh_segments()

        InfoBar.success("Saved", f"Wrote {out_path} + log", parent=self)

    def _check_breakpoints(self) -> None:
        if not self.segments:
            InfoBar.info("No segments", "Add at least one segment.", parent=self)
            return
        ctx = self.context_spin.value()
        enriched = self._materialize_all(ctx)
        if enriched is None:
            return
        # Use enriched data in the UI preview just like post-export.
        self.segments = enriched
        self._refresh_segments()
        InfoBar.success("Checked", "Breakpoint context ready (no files written).", parent=self)

    def _materialize_all(self, ctx: int) -> Optional[List[Dict[str, object]]]:
        used_sources = {seg["src"] for seg in self.segments}
        for src in used_sources:
            path = self._get_path_for_source(src)  # type: ignore[arg-type]
            if path is None or not path.exists():
                InfoBar.error(
                    "Missing file",
                    f"Set {'target' if src == 't' else 'query'} FASTA path before continuing.",
                    parent=self,
                )
                return None
        for src in used_sources:
            if not self.indexes.get(src):
                self._load_index_for_source(src)
                if not self.indexes.get(src):
                    return None
        handles: Dict[Path, object] = {}
        enriched: List[Dict[str, object]] = []
        try:
            for seg in self.segments:
                enriched.append(self._materialize_segment(seg, ctx, handles))
        except Exception as exc:  # pragma: no cover
            InfoBar.error("Operation failed", str(exc), parent=self)
            return None
        finally:
            for fh in handles.values():
                try:
                    fh.close()
                except Exception:
                    pass
        return enriched

    def _build_log_lines(
        self,
        enriched: List[Dict[str, object]],
        merged_len: int,
        out_path: Path,
        suggested_name: str,
        target_name: str,
        query_name: str,
        ctx: int,
        as_html: bool = False,
    ) -> List[str]:
        wrap = _html_wrap if as_html else _md_wrap
        diff_fn = _highlight_diff_html if as_html else _highlight_diff_md
        junction_fn = _junction_preview_html if as_html else _junction_preview_md
        seg_fmt = (
            lambda seg: html.escape(self._format_segment_label(seg))
            if as_html
            else self._format_segment_label(seg)
        )
        rc_count = sum(1 for seg in enriched if seg.get("reverse"))

        if as_html:
            header = [
                "# Manual stitch log",
                "",
                f"- Target FASTA: {html.escape(self.target_path.text() or 'N/A')}",
                f"- Query FASTA: {html.escape(self.query_path.text() or 'N/A')}",
                f"- Target sequence: {html.escape(target_name)}",
                f"- Query sequence: {html.escape(query_name)}",
                f"- Reverse-complemented segments: {rc_count}",
                f"- Output FASTA: {html.escape(str(out_path)) if out_path != Path('N/A') else 'N/A'}",
                f"- Output name: {html.escape(suggested_name)}",
                f"- Segments: {len(enriched)}",
                f"- Merged length: {merged_len:,} bp",
                "",
                "## Segments",
            ]
            seg_line_tpl = "- [{idx}] {src_label}({src}) {name}{rev} {start}-{end} length {length:,} bp"
            junction_title = f"## Breakpoint contexts (each {ctx}bp)"
            preview_label = "  - Junction preview: {preview}"
            prior_right_label = "  - Prior right flank ({src}:{name} {s}-{e}): {seq}"
            next_left_label = "  - Next left flank ({src}:{name} {s}-{e}): {seq}"
            next_left_pre_label = "  - Next left pre-flank ({src}:{name} {s}-{e}): {seq}"
            left_match_label = "  - Left flanks match (last {len}bp): {seq}"
            left_diff_label = "  - Left flanks differ (last {len}bp, diff highlighted):"
            prior_right_diff_label = "    Prior right flank ({src}:{name} {s}-{e}): {seq}"
            next_left_diff_label = "    Next left pre-flank ({src}:{name} {s}-{e}): {seq}"
            no_left_label = "  - No comparable left context."
            right_match_label = "  - Right flanks match (first {len}bp): {seq}"
            right_diff_label = "  - Right flanks differ (first {len}bp, diff highlighted):"
            prior_right_r_label = "    Prior right flank ({src}:{name} {s}-{e}): {seq}"
            next_left_r_label = "    Next left flank ({src}:{name} {s}-{e}): {seq}"
            no_right_label = "  - No comparable right context."
            seg_label = lambda src: "Target" if src == "t" else "Query"
        else:
            header = [
                "# 手动拼接记录",
                "",
                f"- 目标 FASTA: {self.target_path.text() or 'N/A'}",
                f"- 查询 FASTA: {self.query_path.text() or 'N/A'}",
                f"- 目标序列: {target_name}",
                f"- 查询序列: {query_name}",
                f"- 反向互补: {rc_count} 个片段勾选",
                f"- 输出 FASTA: {out_path}",
                f"- 输出序列名: {suggested_name}",
                f"- 段数: {len(enriched)}",
                f"- 合并长度: {merged_len:,} bp",
                "",
                "## 选取的片段",
            ]
            seg_line_tpl = "- [{idx}] {src_label}({src}) {name}{rev} {start}-{end} 长度 {length:,} bp"
            junction_title = f"## 断点上下文（各取 {ctx}bp）"
            preview_label = "  - 断点上下文: {preview}"
            prior_right_label = "  - 前段右侧 ({src}:{name} {s}-{e}): {seq}"
            next_left_label = "  - 后段左侧 ({src}:{name} {s}-{e}): {seq}"
            next_left_pre_label = "  - 后段左侧前序 ({src}:{name} {s}-{e}): {seq}"
            left_match_label = "  - 左侧一致（末{len}bp）: {seq}"
            left_diff_label = "  - 左侧不一致（末{len}bp，高亮为差异位点）:"
            prior_right_diff_label = "    前段右侧 ({src}:{name} {s}-{e}): {seq}"
            next_left_diff_label = "    后段左侧前序 ({src}:{name} {s}-{e}): {seq}"
            no_left_label = "  - 左侧无可比对的上下文。"
            right_match_label = "  - 右侧一致（首{len}bp）: {seq}"
            right_diff_label = "  - 右侧不一致（首{len}bp，高亮为差异位点）:"
            prior_right_r_label = "    前段右侧 ({src}:{name} {s}-{e}): {seq}"
            next_left_r_label = "    后段左侧 ({src}:{name} {s}-{e}): {seq}"
            no_right_label = "  - 右侧无可比对的上下文。"
            seg_label = lambda src: "目标" if src == "t" else "查询"

        log_lines = list(header)
        for idx, seg in enumerate(enriched):
            log_lines.append(
                seg_line_tpl.format(
                    idx=idx,
                    src_label=seg_label(seg["src"]),
                    src=seg["src"],
                    name=html.escape(str(seg["name"])) if as_html else seg["name"],
                    rev=" (RC)" if seg.get("reverse") else "",
                    start=seg["start"],
                    end=seg["end"],
                    length=len(seg["seq"]),
                )
            )

        if len(enriched) >= 2:
            log_lines.append("")
            log_lines.append(junction_title)
            breakpoint_pos = []
            offset = 0
            for i in range(len(enriched) - 1):
                left = enriched[i]
                right = enriched[i + 1]
                offset += len(left["seq"])
                breakpoint_pos.append((i, offset))
                preview = junction_fn(left["right_before"], right["left_after"], ctx)
                log_lines.append(
                    f"- [{i}] {seg_fmt(left)} -> {seg_fmt(right)}"
                )
                log_lines.append(preview_label.format(preview=preview))

                l_rb_len = len(left["right_before"])
                r_la_len = len(right["left_after"])
                r_lb_len = len(right["left_before"])
                l_ctx_start = max(0, left["end"] - l_rb_len)
                l_ctx_end = left["end"]
                r_ctx_start = right["start"]
                r_ctx_end = right["start"] + r_la_len
                log_lines.append(
                    prior_right_label.format(
                        src=left["src"], name=left["name"], s=l_ctx_start, e=l_ctx_end, seq=left["right_before"]
                    )
                )
                log_lines.append(
                    next_left_label.format(
                        src=right["src"], name=right["name"], s=r_ctx_start, e=r_ctx_end, seq=right["left_after"]
                    )
                )
                if r_lb_len:
                    rb_start = max(0, right["start"] - r_lb_len)
                    rb_end = right["start"]
                    log_lines.append(
                        next_left_pre_label.format(
                            src=right["src"], name=right["name"], s=rb_start, e=rb_end, seq=right["left_before"]
                        )
                    )

                left_len = min(50, len(left["right_before"]), len(right["left_before"]))
                if left_len:
                    prev_slice = left["right_before"][-left_len:]
                    curr_slice = right["left_before"][-left_len:]
                    rb_start = max(0, right["start"] - len(right["left_before"]))
                    rb_end = right["start"]
                    if prev_slice == curr_slice:
                        log_lines.append(
                            left_match_label.format(len=left_len, seq=wrap(prev_slice, "green"))
                        )
                    else:
                        hi_prev, hi_curr = diff_fn(prev_slice, curr_slice)
                        log_lines.append(
                            left_diff_label.format(len=left_len)
                        )
                        log_lines.append(
                            prior_right_diff_label.format(
                                src=left["src"], name=left["name"], s=l_ctx_start, e=l_ctx_end, seq=hi_prev
                            )
                        )
                        log_lines.append(
                            next_left_diff_label.format(
                                src=right["src"], name=right["name"], s=rb_start, e=rb_end, seq=hi_curr
                            )
                        )
                else:
                    log_lines.append(no_left_label)

                right_len = min(50, len(left["right_after"]), len(right["left_after"]))
                if right_len:
                    prev_slice = left["right_after"][:right_len]
                    curr_slice = right["left_after"][:right_len]
                    la_start = left["end"]
                    la_end = left["end"] + len(left["right_after"])
                    ra_start = right["start"]
                    ra_end = right["start"] + len(right["left_after"])
                    if prev_slice == curr_slice:
                        log_lines.append(
                            right_match_label.format(len=right_len, seq=wrap(prev_slice, "green"))
                        )
                    else:
                        hi_prev, hi_curr = diff_fn(prev_slice, curr_slice)
                        log_lines.append(
                            right_diff_label.format(len=right_len)
                        )
                        log_lines.append(
                            prior_right_r_label.format(
                                src=left["src"], name=left["name"], s=la_start, e=la_end, seq=hi_prev
                            )
                        )
                        log_lines.append(
                            next_left_r_label.format(
                                src=right["src"], name=right["name"], s=ra_start, e=ra_end, seq=hi_curr
                            )
                        )
                else:
                    log_lines.append(no_right_label)
            log_lines.append("")
            bp_title = "## Breakpoint positions in merged sequence" if as_html else "## 断点在合并序列中的位置"
            log_lines.append(bp_title)
            for i, pos in breakpoint_pos:
                line = f"- breakpoint [{i}] at {pos:,} (0-based, after segment {i})"
                if not as_html:
                    line = f"- 断点 [{i}] 在 {pos:,}（0基，位于第 {i} 段之后）"
                log_lines.append(line)
        if as_html:
            safe_lines: List[str] = []
            for line in log_lines:
                if "<span" in line:
                    safe_lines.append(line)
                else:
                    safe_lines.append(html.escape(line))
            return safe_lines
        return log_lines


# -----------------------------
# Main window
# -----------------------------


class GapNeedleWindow(FluentWindow):
    def __init__(self):
        super().__init__()
        base_font = apply_theme(self)
        self.state = AppState()
        self.setWindowTitle("GapNeedle · Fluent GUI")
        self.resize(WINDOW_WIDTH, WINDOW_HEIGHT)
        try:
            self.navigationInterface.setFont(base_font)
        except Exception:
            pass

        self.align_page = self._wrap_page(AlignPage(self.state))
        self.manual_page = self._wrap_page(ManualStitchPage(self.state))
        self.viewer_widget = AlignmentViewer()
        self.viewer_page = self._wrap_page(self.viewer_widget)
        self._apply_font_and_size(self.align_page.widget(), base_font)
        self._apply_font_and_size(self.manual_page.widget(), base_font)
        self._apply_font_and_size(self.viewer_page.widget(), base_font)

        # Register pages with navigation
        self.addSubInterface(
            self.align_page,
            FluentIcon.SEND,
            "Align",
            position=NavigationItemPosition.TOP,
        )

        item_viewer = self.addSubInterface(
            self.viewer_page,
            FluentIcon.VIEW,
            "PAF viewer",
            NavigationItemPosition.TOP,
        )
        item_manual = self.addSubInterface(
            self.manual_page,
            FluentIcon.EDIT,
            "Manual stitch",
            NavigationItemPosition.TOP,
        )
        self.stackedWidget.setCurrentWidget(self.align_page)
        self._tune_navigation(base_font)
        self._style_nav_item(item_manual, base_font)
        self._style_nav_item(item_viewer, base_font)

    def _tune_navigation(self, base_font: QtGui.QFont) -> None:
        nav_font = QtGui.QFont(base_font)
        nav_font.setPointSize(NAV_FONT_SIZE)
        try:
            self.navigationInterface.panel.setMinimumExpandWidth(NAV_EXPAND_WIDTH)
            self.navigationInterface.panel.setExpandWidth(NAV_EXPAND_WIDTH)
            self.navigationInterface.setMinimumWidth(NAV_MIN_WIDTH)
        except Exception:
            pass
        for btn in self.navigationInterface.findChildren(NavigationPushButton):
            btn.setFont(nav_font)
            btn.setMinimumHeight(NAV_BUTTON_HEIGHT)
            btn.setMinimumWidth(NAV_EXPAND_WIDTH - 20)
            try:
                icon_dim = max(NAV_ICON_SIZE, min(ICON_SIZE, max(16, NAV_BUTTON_HEIGHT - 8)))
                btn.setIconSize(QtCore.QSize(icon_dim, icon_dim))
            except Exception:
                pass

    def _style_nav_item(self, item, base_font: QtGui.QFont) -> None:
        try:
            widget = item.itemWidget
        except Exception:
            return
        nav_font = QtGui.QFont(base_font)
        nav_font.setPointSize(NAV_FONT_SIZE)
        widget.setFont(nav_font)
        widget.setMinimumHeight(NAV_BUTTON_HEIGHT)
        widget.setMaximumHeight(NAV_BUTTON_HEIGHT)
        widget.setMinimumWidth(NAV_EXPAND_WIDTH - 20)
        widget.setMaximumWidth(NAV_EXPAND_WIDTH - 12)

    def _apply_font_and_size(self, widget: QtWidgets.QWidget, font: QtGui.QFont) -> None:
        """Ensure controls inherit scaled font and reasonable heights."""
        if widget is None:
            return
        for child in widget.findChildren(QtWidgets.QWidget):
            child.setFont(font)
            if isinstance(child, (PrimaryPushButton, PushButton, QtWidgets.QAbstractButton)):
                child.setMinimumHeight(CONTROL_HEIGHT)
                if hasattr(child, "setIconSize"):
                    icon_dim = min(ICON_SIZE, max(16, CONTROL_HEIGHT - 8))
                    child.setIconSize(QtCore.QSize(icon_dim, icon_dim))
            elif isinstance(
                child, (LineEdit, ComboBox, SpinBox, QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)
            ):
                child.setMinimumHeight(CONTROL_HEIGHT)

    def _wrap_page(self, widget: QtWidgets.QWidget) -> ScrollArea:
        """Wrap page with a scroll area so full-screen/resize keeps layout intact."""
        area = ScrollArea(self)
        area.setWidgetResizable(True)
        area.setWidget(widget)
        area.setObjectName(widget.objectName() or "page")
        return area


def apply_theme(window: FluentWindow) -> QtGui.QFont:
    setTheme(Theme.LIGHT)
    setThemeColor(QtGui.QColor("#0f766e"))
    # Override global font
    available = {f.lower() for f in QtGui.QFontDatabase().families()}
    chosen = next((c for c in FONT_CANDIDATES if c.lower() in available), "Segoe UI")
    base_font = QtGui.QFont(chosen, BASE_FONT_SIZE)
    app = QtWidgets.QApplication.instance()
    if app:
        app.setFont(base_font)
    return base_font


def launch() -> None:
    import sys

    app = QtWidgets.QApplication.instance()
    if app is None:
        QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
        QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
        app = QtWidgets.QApplication(sys.argv)
    window = GapNeedleWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":  # pragma: no cover
    launch()
