"""
GapFillet 图形界面入口。

提供 `launch()` 便捷函数，启动基于 Tkinter 的双标签页界面：
- 对齐：选择 FASTA 与序列名（带翻页/搜索），设置线程数、preset、反向互补并启动 mappy 比对。
- 手动拼接：沿用上一页的选择上下文，交互式添加片段坐标、预览断点并导出合并结果。
"""

from .app import launch

__all__ = ["launch"]
