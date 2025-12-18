from pathlib import Path
from typing import Dict, Optional, Set

try:
    import mappy as mp
except ImportError as exc:
    raise ImportError(
        "mappy (minimap2 Python bindings) is required. Install with `pip install mappy`."
    ) from exc


def read_fasta_sequences(path: Path, select_names: Optional[Set[str]] = None) -> Dict[str, str]:
    """
    Load sequences from FASTA into a dict, using mappy.fastx_read for speed.
    If select_names is provided, stop early once all requested names are loaded to avoid reading the entire file.
    """
    path = Path(path)
    sequences: Dict[str, str] = {}
    for name, seq, _ in mp.fastx_read(str(path), read_comment=False):
        key = name.split()[0]
        if select_names is None or key in select_names:
            sequences[key] = seq
            if select_names and len(sequences) == len(select_names):
                break
    return sequences


def write_fasta(records: Dict[str, str], path: Path) -> None:
    """Write a simple FASTA file from a name->sequence dictionary."""
    path = Path(path)
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]
