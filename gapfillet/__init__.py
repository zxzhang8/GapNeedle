"""
GapFillet: utilities for gap filling and stitching genome assemblies with minimap2.

Key entrypoints
---------------
- run_minimap2_alignment: thin wrapper around minimap2 with safe output handling.
- suggest_overlaps: read PAF and summarize candidate joins.
- stitch_from_paf: stitch two sequences using a chosen alignment from PAF.
- GapFillet: convenience facade combining the above for interactive use.
"""

from .aligner import run_minimap2_alignment
from .stitcher import (
    GapFillet,
    manual_stitch_by_coordinates,
    suggest_overlaps,
    stitch_from_paf,
    telomere_presence,
)
from .gui import launch
from .io import read_fasta_sequences, write_fasta

__all__ = [
    "GapFillet",
    "run_minimap2_alignment",
    "suggest_overlaps",
    "stitch_from_paf",
    "telomere_presence",
    "manual_stitch_by_coordinates",
    "read_fasta_sequences",
    "write_fasta",
    "launch",
]
