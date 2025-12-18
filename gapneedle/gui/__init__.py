"""
Entry point for the GapNeedle GUI.

Exposes a convenient `launch()` function that starts the two-tab Tkinter interface:
- Align: pick FASTA files and sequence names (with paging/search), set threads/preset/reverse-complement, then run mappy.
- Manual stitch: reuse selections from the Align tab, interactively add segment coordinates, preview junctions, and export the merged result.
"""

from .app import launch

__all__ = ["launch"]
