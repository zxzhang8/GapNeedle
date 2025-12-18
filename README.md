GapNeedle
=========

一个面向基因组拼接 gap filling 的轻量 Python 库，封装了 minimap2 对单条序列的比对、PAF 解析和简易拼接逻辑，兼顾命令行和 Jupyter 使用体验。

安装
----
- 依赖：Python 3.9+，`mappy`（已安装）。
- [mappy](https://github.com/lh3/minimap2) 安装

  GapNeedle的比对依赖minimap2，mappy是minimap2的Python bindings to a sub set of C APIs。其可以通过 `pip install mappy` 或 `conda install -c bioconda mappy`安装。
- 本地开发/临时使用：在项目根目录执行
  ```bash
  pip install -e .
  ```
  之后可在任意目录使用 `python -m gapfillet.cli ...` 或 `gapfillet ...`。

主要功能
--------
- `run_minimap2_alignment`: 调用 minimap2 输出 PAF，参数与官方 CLI 保持一致（`-x`, `-t` 等），默认只抽取指定序列避免全基因组比对。
- `suggest_overlaps`: 解析 PAF，列出候选重叠片段（按 overlap 长度排序）。
- `stitch_from_paf`: 读取 PAF 挑选一个 alignment，将查询序列的左右延伸与目标序列拼接，strand 自动处理，结果输出 FASTA。
- `GapFillet.scan_gaps`: 快速扫描 FASTA 中长度>=N 的 `N` 区域，辅助 gap 填补决策。
- `telomere_presence` / `GapFillet.check_telomere`: 检查指定序列两端给定窗口内是否含连续端粒基序（默认 CCCTAA 或互补，至少连续15次）。
- `manual_stitch_by_coordinates`: 手动输入来源序列与坐标逐段拼接，支持打印侧翼序列和拼接边界一致性提示。

命名与输出
----------
- 默认输出目录/文件名：`{query}_vs_{target}/{query}_vs_{target}.{preset}.paf`，不同 preset 会写入不同文件；如存在则自动跳过并提示。
- 拼接结果默认写入 `{paf}.stitched.fasta`，序列名为 `{target}+{query}`。

流程示例（命令 + 说明 + 例子）
------------------------------

1) 生成 PAF：将查询序列比对到目标序列  
   - 输入：两个 FASTA、各自的序列名。  
   - 行为：调用 minimap2（mappy），默认只抽取指定单条序列；输出 `{query}_vs_{target}/{query}_vs_{target}.{preset}.paf`（preset 自动写入文件名），如果同名文件存在则直接跳过并提示。  
   - CLI 例子：
     ```bash
     python -m gapneedle.cli align ref.fa qry.fa chr1 chr1 --threads 8 --preset asm10
     # 若已有 chr1_vs_chr1/chr1_vs_chr1.asm10.paf 则跳过运行
     ```
   - Python / Jupyter 例子：
     ```python
     from gapneedle import GapFillet
     gf = GapFillet()
     run = gf.align(
         target_fasta="ref.fa",
         query_fasta="qry.fa",
         target_seq="chr1",
         query_seq="chr1",
         threads=8,
         preset="asm10",
     )
     print(run.output_path)  # 默认 chr1_vs_chr1/chr1_vs_chr1.asm10.paf
     ```

2) 拼接：读取 PAF，选择重叠，输出拼接结果  
   - 输入：PAF、两个 FASTA、序列名。自动检查 PAF；不存在时可通过 `auto_align=True` 先跑对齐。  
   - 行为：列出候选重叠（overlap/strand/坐标/identity），可交互选择编号并确认；自动处理正反向链。拼接后写入 `{paf}.stitched.fasta`。若左/右侧有额外延伸会附加到目标序列两端。  
   - CLI 例子：
     ```bash
     # 交互选择候选（默认 selection 留空即交互）
     python -m gapneedle.cli stitch chr1_vs_chr1/chr1_vs_chr1.asm10.paf ref.fa qry.fa chr1 chr1
     # 或指定候选编号并输出到自定义路径
     python -m gapneedle.cli stitch chr1_vs_chr1/chr1_vs_chr1.asm10.paf ref.fa qry.fa chr1 chr1 \
       --selection 0 --output chr1_merge.fa
     ```
   - Python / Jupyter 例子：
     ```python
     merged = gf.stitch(
         run=run,                       # 直接复用 align 返回的对象，可省去 fasta/序列名/PAF 路径
         selection=None,                # None=交互选择；给整数则自动选
         interactive=True,              # Notebook 下可输入编号
     )
     print(merged)  # 结果 FASTA 路径
     ```

3) 其它辅助：gap 扫描  
   - 作用：快速列出 FASTA 中长度 ≥ N 的 `N` 区域，辅助决定填补位置。  
   - CLI：`python -m gapfillet.cli scan-gaps assembly.fa --min-gap 20`  
   - Python：
     ```python
     gaps = gf.scan_gaps("assembly.fa", min_gap=20)
     # gaps: [(seq_name, start, end), ...]，坐标为 0-based 半开区间
     ```

4) 端粒基序检测  
   - 作用：判断序列两端给定窗口内是否含连续端粒基序（默认 CCCTAA，自动考虑反向互补，至少连续15次才算命中）。  
   - Python：
     ```python
     from gapneedle import telomere_presence
     left, right = telomere_presence("assembly.fa", "chr1", window=10_000_000, min_repeats=15)
     # 或使用封装在 GapNeedle 中的便捷方法
     gf.check_telomere("assembly.fa", "chr1", window=5_000_000, motif="CCCTAA", min_repeats=15)
     # 若需要细节（连续匹配位置/重复数/覆盖长度/变异碱基数），使用 telomere_details
     from gapneedle.stitcher import telomere_details
     left_info, right_info = telomere_details("assembly.fa", "chr1", window=5_000_000, min_repeats=15)
     # info 示例: {'has': True, 'matches': [(0,6), ...], 'span_len': 120, 'mutated_bases': 0, 'repeat_count': 20}
     ```

5) 手动坐标拼接  
   - 作用：在交互模式下按坐标逐段选取目标/查询序列并拼接，可打印两侧 ~200bp 序列并检查相邻段边界是否一致，结束后汇总显示每个断点两侧序列（带颜色高亮和断点符号），最后询问输出文件/序列名并保存，同时生成同前缀的日志文件（.md）记录所有片段坐标。  
   - Python：
     ```python
     from gapneedle import manual_stitch_by_coordinates
     merged = manual_stitch_by_coordinates(
         "target.fa",
         "query.fa",
         "chr1",
         "chr1",
         context=200,                    # 断点上下文长度，同时用于汇总
     )
     # 按提示输入 t/q 与 start-end，依次打印断点左右两侧序列（先左断点，确认后给出右断点），结束后汇总显示各断点，再输入保存文件路径和序列名；输入 stop_token(默认 x，不能与 t/q 重复) 结束。保存时会同时写出同前缀的 .md 日志，记录每段的来源和坐标。
     ```

图形界面（GUI）
---------------
- 双标签页：`对齐`（选择 FASTA，分页搜索序列名，设置线程/preset/反向互补后启动 mappy）与 `手动拼接`（继承上一页选择，输入坐标即时预览断点一致性并导出）。
- 设计：暗色卡片 + 荧光绿强调，日志滚动面板实时输出；长序列列表支持搜索与翻页，避免一次性加载过多行。
- 启动：
  ```bash
  # 需要额外安装 PyQt5 与 PyQt-Fluent-Widgets（未写入核心依赖）
  pip install PyQt5 PyQt-Fluent-Widgets

  python -m gapneedle.gui.app            # 或 Python 中：
  # from gapneedle.gui import launch; launch()
  ```
- 字体/尺寸自定义：修改 `gapfillet/gui/theme.py` 中的 FONT_CANDIDATES、BASE_FONT_SIZE、WIDGET_FONT_SIZE 等即可统一调整字体和按钮大小。

Jupyter 友好
------------
- 所有高层 API 接受 `say` 回调，自定义消息输出（默认 `print`）。
- `stitch_from_paf` 支持 `interactive=True` 时在 notebook 中通过输入编号选择候选拼接。

注意事项
--------
- 需要本地安装 minimap2 的 Python 绑定 `mappy`（官方提供，`pip install mappy`），内部通过 mappy 调用 minimap2 核心算法。
- 拼接逻辑为快速合并：保留目标序列主体，在左右两端用查询序列的额外覆盖延伸。若需要更复杂的共识或拼接策略，可在此基础上扩展。
