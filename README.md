GapNeedle (C++/Qt)
==================

English README · [中文文档](README_cn.md)

<img src="img/logo.svg" alt="GapNeedle Logo" width="180" />

GapNeedle is a desktop-first gap-filling tool implemented in C++17 with Qt6 Widgets.
It focuses on efficient single-sequence workflows: **align -> inspect -> stitch**.

Current core capabilities:
- FASTA IO (read/write), reverse-complement, and indexed slice access via `.fai` (auto-build on demand).
- PAF parsing/filtering and overlap suggestion.
- Query index -> target index mapping from `cg:Z` CIGAR.
- Segment-based stitch service, gap scanning, and telomere motif check.

Contents
--------
- Quickstart
- Build
- Run
- GUI Modules
- CLI
- Outputs
- Notes
- Current Limits

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
  - Select target/query FASTA and sequence names (names loaded in background).
  - Configure preset/threads/reverse-complement flags.
  - Uses cache path under app-relative `cache/alignments/...` and reuses existing PAF when available.
  - Syncs alignment context to PAF Viewer and Manual Stitch after success.
- **PAF Viewer**
  - Load/filter/sort records (`mapQ` filter and multiple sort keys).
  - Highlight overlaps on query intervals.
  - Map query index to target index using `cg:Z` CIGAR and show detailed mapping diagnostics.
  - Mapping is unavailable when `cg:Z` is missing in the selected record.
- **Manual Stitch**
  - Build segment list from target/query/extra FASTA sources.
  - Edit segment order (add/remove/move/resume) and optional reverse-complement per segment.
  - Run background breakpoint checks with flank-diff visualization.
  - Export merged FASTA and JSON session log.
  - Load session from JSON, with fallback support for legacy markdown logs.
- **FASTA Search**
  - Search subsequences in FASTA and list hit coordinates.
  - Matching is exact substring matching (no mismatch-tolerant search yet).

CLI
---
`gapneedle_cli` supports:
- `align`
  - Required: `--target-fasta --query-fasta --target-seq --query-seq`
  - Optional: `--output --preset --threads`
- `stitch`
  - Required: `--target-fasta --query-fasta --output`
  - Repeatable segment: `--segment src:name:start:end[:rc]`
  - Optional: `--output-name`
- `scan-gaps`
  - Required: `--target-fasta`
  - Optional: `--min-gap`
- `check-telomere`
  - Required: `--target-fasta --seq-name`

Use `--cmd <name>` with corresponding options.

Outputs
-------
- Alignment PAF:
  - GUI align flow: app-relative `cache/alignments/...`
  - Core aligner default path (when output is omitted): `resources/<query_vs_target>/<...>.paf`
- Stitch output FASTA: user-selected path
- Stitch session log: `<output>.session.json`
- Viewer/search modules: no fixed output artifact unless exported/copied by user actions.

Notes
-----
- `GAPNEEDLE_USE_MINIMAP2=ON` expects minimap2 sources under `third_party/minimap2/`.
- When minimap2 sources are available, source-based minimap2 bridge is built and used.
- When minimap2 is unavailable (or disabled), build falls back to a stub aligner:
  - If requested output PAF already exists and reuse is enabled, align can return cached result.
  - Otherwise align fails with a minimap2 integration error.
- Query->target coordinate mapping depends on `cg:Z` in PAF records.

Current Limits
--------------
- Workflow is still manual-decision first for stitching (no automatic stitch planner).
- FASTA Search is exact-match only.
- Automated test coverage is currently lightweight.
