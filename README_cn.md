GapNeedle（C++/Qt）
===================

中文文档 · [English README](README.md)

GapNeedle 是一个基于 C++17 + Qt6 Widgets 的桌面 gap-filling 工具，
核心流程是 **对齐 -> 查看 -> 拼接**。

当前核心能力：
- FASTA 读写、反向互补、基于 `.fai` 的索引切片读取（按需自动建索引）。
- PAF 解析/过滤与重叠候选建议。
- 基于 `cg:Z` CIGAR 的 query 坐标 -> target 坐标映射。
- 基于片段的拼接服务、gap 扫描与端粒 motif 检查。

目录
----
- 快速开始
- 构建
- 运行
- GUI 模块
- CLI
- 输出
- 说明
- 当前边界

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
  - 选择 target/query FASTA 与序列名（序列名后台加载）
  - 设置 preset/线程/反向互补
  - 使用应用相对目录 `cache/alignments/...` 进行缓存，并在可用时复用已有 PAF
  - 对齐成功后自动同步上下文到 PAF Viewer 与 Manual Stitch
- **PAF Viewer**
  - 记录过滤、排序（含 `mapQ` 过滤与多字段排序）
  - query 区间重叠标记
  - 基于 `cg:Z` CIGAR 执行 query->target 坐标映射并展示详细诊断信息
  - 若记录缺失 `cg:Z`，则无法进行精确坐标映射
- **Manual Stitch**
  - 从 target/query/extra FASTA 组装片段
  - 片段编辑（添加/删除/上下移动/回填编辑），支持片段级反向互补
  - 后台执行断点检查并可视化 flank 差异
  - 导出合并 FASTA 与 JSON 会话日志
  - 支持从 JSON 会话加载，并兼容部分 legacy markdown 会话日志
- **FASTA Search**
  - FASTA 子序列检索与坐标结果展示
  - 当前为精确子串匹配（暂不支持错配容忍搜索）

CLI
---
`gapneedle_cli` 支持：
- `align`
  - 必需：`--target-fasta --query-fasta --target-seq --query-seq`
  - 可选：`--output --preset --threads`
- `stitch`
  - 必需：`--target-fasta --query-fasta --output`
  - 可重复片段参数：`--segment src:name:start:end[:rc]`
  - 可选：`--output-name`
- `scan-gaps`
  - 必需：`--target-fasta`
  - 可选：`--min-gap`
- `check-telomere`
  - 必需：`--target-fasta --seq-name`

通过 `--cmd <命令>` 调用。

输出
----
- 对齐 PAF：
  - GUI 对齐流程：应用相对目录 `cache/alignments/...`
  - Core 对齐器默认输出（未显式给 `--output` 时）：`resources/<query_vs_target>/<...>.paf`
- 拼接输出 FASTA：用户选择路径
- 拼接会话日志：`<output>.session.json`
- Viewer/Search 模块默认不生成固定产物（除用户显式导出/复制行为）。

说明
----
- 当 `GAPNEEDLE_USE_MINIMAP2=ON` 时，要求 `third_party/minimap2/` 下存在 minimap2 源码。
- 当 minimap2 源码可用时，会构建并使用源码桥接层执行真实对齐。
- 当 minimap2 不可用（或关闭）时，会回退到桩实现：
  - 若请求输出路径已有 PAF 且允许复用，可直接返回缓存结果。
  - 否则 `align` 会报 minimap2 集成不可用错误。
- query->target 坐标映射依赖 PAF 记录中的 `cg:Z` 字段。

当前边界
--------
- 拼接流程仍以人工决策为主（暂无自动拼接规划器）。
- FASTA Search 当前仅支持精确匹配。
- 自动化测试覆盖目前较轻量。
