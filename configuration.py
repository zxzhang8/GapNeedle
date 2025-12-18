"""
User-tunable hyperparameters for GapNeedle.

Place this file at the project root and tweak values as needed; other modules can import it (e.g., GUI theme or default paths).
All paths can be str or pathlib.Path.
"""

from pathlib import Path

# ==== File dialog defaults ====
# Starting directory when opening FASTA files; set to None to use the OS default/last visited dir.
DEFAULT_BROWSE_DIR: str | Path | None = None

# Starting directory when saving stitched results or other outputs; set to None to follow dialog defaults.
DEFAULT_SAVE_DIR: str | Path | None = None

# Default directory for alignment outputs such as PAF files.
ALIGN_OUTPUT_DIR = str(Path(__file__).resolve().parent / "resources")

# ==== GUI theme/layout hyperparameters (mirrors gapneedle/gui/theme.py) ====
# Font candidates ordered by priority; first available is used.
FONT_CANDIDATES = [
    "Microsoft YaHei UI",
    "Microsoft YaHei",
    "PingFang SC",
    "Noto Sans CJK SC",
    "WenQuanYi Zen Hei",
    "SimHei",
    "Arial",
    "Helvetica",
    "TkDefaultFont",
]

# Font size scaling; 1.0 keeps original size, >1 enlarges, <1 shrinks.
FONT_SCALING = 2

# Control size scaling; affects button height, icon size, etc.
UI_SCALING = 1.5

# Initial window size (pixels).
WINDOW_WIDTH = 2880
WINDOW_HEIGHT = 1620

# Font size presets (points), derived from FONT_SCALING.
BASE_FONT_SIZE = int(12 * FONT_SCALING)
ACCENT_FONT_SIZE = int(12 * FONT_SCALING)
HERO_FONT_SIZE = int(16 * FONT_SCALING)
WIDGET_FONT_SIZE = int(16 * FONT_SCALING)

# Control height and icon size (pixels), impacted by UI_SCALING.
CONTROL_HEIGHT = int(36 * UI_SCALING)
ICON_SIZE = int(20 * UI_SCALING)

# Padding (pixels) for buttons, cards, panels, etc.
BUTTON_PADDING = int(12 * UI_SCALING)
CARD_PADDING = int(12 * UI_SCALING)
PANEL_PADDING = int(14 * UI_SCALING)

# List pagination size (items per page).
LIST_PAGE_SIZE = 12

# Navigation bar sizes and font/icon settings.
NAV_MIN_WIDTH = int(72 * UI_SCALING)
NAV_EXPAND_WIDTH = int(220 * UI_SCALING)
NAV_BUTTON_HEIGHT = int(44 * UI_SCALING)
NAV_FONT_SIZE = BASE_FONT_SIZE
NAV_ICON_SIZE = int(32 * UI_SCALING)


# Info panel height (pixels) for the manual stitch UI.
INFO_HEIGHT = 500
