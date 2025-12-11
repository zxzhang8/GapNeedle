"""
Global UI tuning knobs for GapFillet GUI.

Adjust sizes and font candidates here to quickly reskin the app without hunting
through the UI code.
"""

# Font candidate priority (first available wins)
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


# Optional font scaling; 1.0 = default.
FONT_SCALING = 2
# Optional UI scaling factor for control heights/widths/icons; 1.0 = default.
UI_SCALING = 1.5

# Initial window size (width, height)
WINDOW_WIDTH = 2880
WINDOW_HEIGHT = 1620

# Size presets (in points)
BASE_FONT_SIZE = int(12 * FONT_SCALING)
ACCENT_FONT_SIZE = int(12 * FONT_SCALING)
HERO_FONT_SIZE = int(16 * FONT_SCALING)
WIDGET_FONT_SIZE = int(16 * FONT_SCALING)
CONTROL_HEIGHT = int(36 * UI_SCALING)
ICON_SIZE = int(20 * UI_SCALING)

# Padding tweaks (pixels)
BUTTON_PADDING = int(12 * UI_SCALING)
CARD_PADDING = int(12 * UI_SCALING)
PANEL_PADDING = int(14 * UI_SCALING)
LIST_PAGE_SIZE = 12

# Navigation tuning
NAV_MIN_WIDTH = int(72 * UI_SCALING)
NAV_EXPAND_WIDTH = int(220 * UI_SCALING)
NAV_BUTTON_HEIGHT = int(44 * UI_SCALING)
NAV_FONT_SIZE = BASE_FONT_SIZE
NAV_ICON_SIZE = int(32 * UI_SCALING)

# Default working directory for file dialogs (None means OS default)
DEFAULT_BROWSE_DIR = "/home/zxzhang/remote/mouse/KM"  # set to "/path/to/data" if you want dialogs to start there
ALIGN_OUTPUT_DIR = "D:/projects/GapFillet/resources"

# MANUAL STITCH UI
INFO_HEIGHT = 500
