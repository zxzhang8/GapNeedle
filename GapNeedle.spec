# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path

block_cipher = None

spec_path = Path(globals().get("__file__", Path.cwd())).resolve()
project_root = spec_path.parent
app_entry = project_root / "GapNeedle" / "gapneedle" / "gui" / "app.py"

datas = []

config_file = project_root / "configuration.py"
if config_file.exists():
    datas.append((str(config_file), "."))

resources_dir = project_root / "resources"
if resources_dir.exists():
    datas.append((str(resources_dir), "resources"))

analysis = Analysis(
    [str(app_entry)],
    pathex=[str(project_root)],
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "PyQt6",
        "PySide6",
        "PySide2",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(analysis.pure, analysis.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    analysis.scripts,
    analysis.binaries,
    analysis.zipfiles,
    analysis.datas,
    [],
    name="GapNeedle",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
