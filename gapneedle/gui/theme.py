"""
Global UI tuning knobs for GapNeedle GUI.

Values are loaded from the project-root `configuration.py` when available,
otherwise fall back to the defaults defined here. This keeps GUI styling and
paths configurable without editing package files.
"""

from importlib import import_module
from pathlib import Path
from typing import Any, Optional

try:
    _CFG = import_module("configuration")
except Exception:
    _CFG = None


def _cfg_value(name: str, default: Any) -> Any:
    cfg = _CFG
    return getattr(cfg, name, default) if cfg else default


# Font candidate priority (first available wins)
FONT_CANDIDATES = _cfg_value(
    "FONT_CANDIDATES",
    [
        "Microsoft YaHei UI",
        "Microsoft YaHei",
        "PingFang SC",
        "Noto Sans CJK SC",
        "WenQuanYi Zen Hei",
        "SimHei",
        "Arial",
        "Helvetica",
        "TkDefaultFont",
    ],
)

# Optional font scaling; 1.0 = default.
FONT_SCALING: float = _cfg_value("FONT_SCALING", 2)
# Optional UI scaling factor for control heights/widths/icons; 1.0 = default.
UI_SCALING: float = _cfg_value("UI_SCALING", 1.5)

# Initial window size (width, height)
WINDOW_WIDTH: int = _cfg_value("WINDOW_WIDTH", 2880)
WINDOW_HEIGHT: int = _cfg_value("WINDOW_HEIGHT", 1620)

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
LIST_PAGE_SIZE: int = _cfg_value("LIST_PAGE_SIZE", 12)

# Navigation tuning
NAV_MIN_WIDTH = int(72 * UI_SCALING)
NAV_EXPAND_WIDTH = int(220 * UI_SCALING)
NAV_BUTTON_HEIGHT = int(44 * UI_SCALING)
NAV_FONT_SIZE = BASE_FONT_SIZE
NAV_ICON_SIZE = int(32 * UI_SCALING)

# Default working directories (None means OS default)
DEFAULT_BROWSE_DIR: Optional[str | Path] = _cfg_value("DEFAULT_BROWSE_DIR", None)
DEFAULT_SAVE_DIR: Optional[str | Path] = _cfg_value("DEFAULT_SAVE_DIR", None)
ALIGN_OUTPUT_DIR = _cfg_value("ALIGN_OUTPUT_DIR", str(Path(__file__).resolve().parent.parent / "resources"))

# MANUAL STITCH UI
INFO_HEIGHT: int = _cfg_value("INFO_HEIGHT", 500)
