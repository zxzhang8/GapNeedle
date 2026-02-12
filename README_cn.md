GapNeedle（C++/Qt）
===================

中文文档 · [English README](README.md)

GapNeedle 是一个基于 C++17 + Qt6 Widgets 的桌面 gap-filling 工具，
核心流程是 **对齐 -> 查看 -> 拼接**。

目录
----
- 快速开始
- 构建
- 运行
- GUI 模块
- CLI
- 输出
- 说明

快速开始
--------
```bash
cmake -S . -B build_cpp -DGAPNEEDLE_BUILD_GUI=ON -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp -j
```

构建
----
依赖要求：
- CMake 3.21+
- 支持 C++17 的编译器
- Qt6（Core/Gui/Widgets）

常用选项：
- `GAPNEEDLE_BUILD_GUI=ON|OFF`
- `GAPNEEDLE_BUILD_CLI=ON|OFF`
- `GAPNEEDLE_BUILD_TESTS=ON|OFF`
- `GAPNEEDLE_USE_MINIMAP2=ON|OFF`

示例：
```bash
# core + cli + tests
cmake -S . -B build_cpp -DGAPNEEDLE_BUILD_GUI=OFF -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp -j
ctest --test-dir build_cpp --output-on-failure

# gui
cmake -S . -B build_cpp_gui -DGAPNEEDLE_BUILD_GUI=ON -DGAPNEEDLE_USE_MINIMAP2=OFF
cmake --build build_cpp_gui -j
```

运行
----
CLI：
```bash
./build_cpp/gapneedle_cli --cmd scan-gaps --target-fasta /path/to/ref.fa --min-gap 10
```

GUI：
```bash
./build_cpp_gui/gapneedle_gui
```

GUI 模块
--------
- **Align**
  - 选择 target/query FASTA 与序列名
  - 设置 preset/线程/反向互补
  - 自动使用应用相对目录下的 `cache/alignments/...` 进行缓存
- **PAF Viewer**
  - 记录过滤、排序、重叠标记
  - 基于 `cg:Z` CIGAR 执行 query->target 坐标映射
- **Manual Stitch**
  - 从 target/query/extra FASTA 组装片段
  - 断点一致性检查
  - 导出合并 FASTA 与 JSON 会话日志
- **FASTA Search**
  - FASTA 子序列检索与坐标结果展示

CLI
---
`gapneedle_cli` 支持：
- `align`
- `stitch`
- `scan-gaps`
- `check-telomere`

通过 `--cmd <命令>` 调用。

输出
----
- 对齐缓存 PAF：应用相对目录 `cache/alignments/...`
- 拼接输出 FASTA：用户选择路径
- 拼接会话日志：`<output>.session.json`

说明
----
- 当 `GAPNEEDLE_USE_MINIMAP2=ON` 时，要求 `third_party/minimap2/` 下存在 minimap2 源码。
- minimap2 桥接代码位于 `src/minimap2_bridge/`，可在此基础上继续完善。
