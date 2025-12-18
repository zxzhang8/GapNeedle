GapNeedle
=========

轻量级 gap filling 工具，封装 minimap2（mappy）比对、PAF 解析和快速拼接逻辑。提供 GUI 与 Python API，专注单序列到单序列的高效对齐、查看与拼接。

目录
----
- 快速上手
- 安装与依赖
- GUI 使用指南
  - 比对（Alignment）
  - 结果查看（Alignment Viewer）
  - 拼接（Stitching）
- 反向互补比对注意事项
- Python API 示例
- 配置与默认输出
- 常见问题

快速上手
--------
```bash
# 一键安装（包含 GUI 依赖）
pip install ".[gui]"

# 运行 GUI
python -m gapneedle.gui.app
```
在 GUI 中完成比对、查看和拼接；或使用 Python API 进行编程调用。

安装与依赖
----------
- Python 3.9+
- 核心依赖：`mappy`（随 `pip install ".[gui]"` 一并安装）
- GUI 依赖：`PyQt5`, `PyQt-Fluent-Widgets`（同上）
若只需核心功能，可安装 `pip install .`；开发模式可用 `pip install -e ".[gui]"`。

GUI 使用指南
-----------
界面包含“对齐”和“手动拼接”两个标签页，日志面板实时输出进度，长序列列表支持分页搜索。

### 比对（Alignment）
- 输入：目标 FASTA、查询 FASTA、各自序列名；可选线程数、preset、是否反向互补查询序列。
- 输出：PAF 文件，默认路径 `{query}.{qseq}_vs_{target}.{tseq}/{...}.paf` 写入 `resources/`（或配置的目录）。

操作示例：
1. 在“对齐”页点击 Browse 选择目标/查询 FASTA，选择序列名。
2. 设置线程（默认 4）、preset（默认 asm20，亦可 asm5/asm10/map-ont/...），勾选是否反向互补查询。
3. 点击“Run alignment”。完成后日志显示 PAF 路径，并缓存选择以便后续拼接。

### 结果查看（Alignment Viewer）
- 比对完成后，点击日志中的 PAF 路径即可在右侧查看记录，或手动加载 PAF。
- 展示：查询/目标长度与坐标、链方向（-> 或 <-）、matches、比对长度、mapq。

### 拼接（Stitching）
- 输入：PAF、目标/查询 FASTA 与序列名；可复用“对齐”页上下文。
- 过程：
  1. 切到“手动拼接”页，点击“Fill from Align tab”自动填充。
  2. 选择 PAF，自动列出候选 overlap（按长度排序）。
  3. 选择候选编号，预览断点上下文（默认各取 200bp），一致绿色，不一致橙色。
  4. 点击 Export，输入输出路径/序列名（默认 `<target>+<query>`），生成 FASTA 和同名前缀的 `.md` 日志，记录片段来源与断点。

反向互补比对注意事项
------------------
- 为什么需要：查询序列方向不确定或已知应为互补链时，先取反向互补可获得正确的比对与后续拼接方向。
- 如何使用：在 “对齐” 页勾选 “Reverse-complement query” 后再运行比对。当前仅支持对查询序列取反向互补，Target 不支持反向互补比对。
- 坐标与一致性：勾选后，Alignment 生成的 PAF 坐标基于反向互补后的查询序列；Stitch 页继续勾选使用时，拼接和断点坐标与 Alignment 结果保持一致（均基于反向互补后的序列）。
- 注意：如果未勾选但真实方向相反，可能覆盖度低或比对失败；若误勾选，坐标会对应互补方向，请确认后再拼接。

Python API 示例
--------------
```python
from gapneedle import GapNeedle

gf = GapNeedle()
run = gf.align("ref.fa", "qry.fa", "chr1", "chr1", threads=8, preset="asm10")
print("PAF:", run.output_path)

merged_path = gf.stitch(run=run, selection=None, interactive=False)  # 自动选第 1 条
print("Merged FASTA:", merged_path)
```

配置与默认输出
-------------
- 配置文件：项目根 `configuration.py`（或包内默认值）。可设置字体、窗口尺寸、默认打开/保存目录、PAF 输出目录等。
- 输出命名：
  - PAF：`{query}.{qseq}_vs_{target}.{tseq}.{preset}.paf`，存放于 `resources/` 或配置目录。
  - 拼接结果：`<paf>.stitched.fasta`（序列名 `<target>+<query>`），或 GUI 指定路径。
  - 拼接日志：与输出 FASTA 同名前缀的 `.md`。

常见问题
--------
- **mappy 导入失败**：先安装 `pip install mappy`（或 `conda install -c bioconda mappy`），确保 minimap2 动态库可用。
- **GUI 配置未生效**：确认已修改 `configuration.py`，并重启 GUI。
- **PAF 未生成**：检查序列名是否存在；若 PAF 已存在且开启复用，会跳过运行（日志提示 “Existing PAF detected”）。
