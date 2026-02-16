GapNeedle (C++/Qt)
==================

English README · [中文文档](README_cn.md)

<img src="img/logo.svg" alt="GapNeedle Logo" width="180" />

GapNeedle is a desktop-first gap-filling tool implemented in C++17 with Qt6 Widgets.
It focuses on efficient single-sequence workflows: **align -> inspect -> stitch**.

Contents
--------
- Quickstart
- Build
- Run
- GUI Modules
- CLI
- Outputs
- Notes

Quickstart
----------
```bash
cmake -S . -B build_cpp -DGAPNEEDLE_BUILD_GUI=ON -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp -j
```

Build
-----
Requirements:
- CMake 3.21+
- C++17 compiler
- Qt6 (Core/Gui/Widgets)

Options:
- `GAPNEEDLE_BUILD_GUI=ON|OFF`
- `GAPNEEDLE_BUILD_CLI=ON|OFF`
- `GAPNEEDLE_BUILD_TESTS=ON|OFF`
- `GAPNEEDLE_USE_MINIMAP2=ON|OFF`

Examples:
```bash
# core + cli + tests
cmake -S . -B build_cpp -DGAPNEEDLE_BUILD_GUI=OFF -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp -j
ctest --test-dir build_cpp --output-on-failure

# gui build
cmake -S . -B build_cpp_gui -DGAPNEEDLE_BUILD_GUI=ON -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp_gui -j
```

Run
---
CLI:
```bash
./build_cpp/gapneedle_cli --cmd scan-gaps --target-fasta /path/to/ref.fa --min-gap 10
```

GUI:
```bash
./build_cpp_gui/gapneedle_gui
```

GUI Modules
-----------
- **Align**
  - Select target/query FASTA and sequence names.
  - Configure preset/threads/reverse-complement flags.
  - Uses cache path under app-relative `cache/alignments/...`.
- **PAF Viewer**
  - Load/filter/sort records.
  - Highlight overlaps.
  - Map query index to target index using `cg:Z` CIGAR.
- **Manual Stitch**
  - Build segment list from target/query/extra FASTA sources.
  - Check breakpoint consistency.
  - Export merged FASTA and JSON session log.
- **FASTA Search**
  - Search subsequences in FASTA and list hit coordinates.

CLI
---
`gapneedle_cli` supports:
- `align`
- `stitch`
- `scan-gaps`
- `check-telomere`

Use `--cmd <name>` with corresponding options.

Outputs
-------
- Alignment cache PAF: app-relative `cache/alignments/...`
- Stitch output FASTA: user-selected path
- Stitch session log: `<output>.session.json`

Notes
-----
- `GAPNEEDLE_USE_MINIMAP2=ON` expects minimap2 sources under `third_party/minimap2/`.
- Current bridge scaffolding exists in `src/minimap2_bridge/`; complete integration as needed.
