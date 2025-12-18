import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

try:
    import mappy as mp
except ImportError as exc:
    raise ImportError(
        "mappy (minimap2 Python bindings) is required. "
        "Install with `pip install mappy` or build from minimap2 source."
    ) from exc

from .io import read_fasta_sequences, reverse_complement, write_fasta


class Minimap2Error(RuntimeError):
    """Raised when minimap2 exits with non-zero status."""


@dataclass
class AlignmentRun:
    target_fasta: Path
    query_fasta: Path
    target_seq: str
    query_seq: str
    preset: str
    output_path: Path
    cmd: List[str]
    skipped: bool


def _safe_part(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)


def default_paf_path(
    target_fasta: Path,
    query_fasta: Path,
    target_seq: str,
    query_seq: str,
    preset: str,
    output_dir: Optional[Path] = None,
    *,
    reverse_query: bool = False,
) -> Path:
    """
    Infer default PAF path. Use FASTA filename + sequence name to avoid collisions when different files share names.
    """
    safe_target = _safe_part(target_seq)
    safe_query = _safe_part(query_seq) + ("_rc" if reverse_query else "")
    safe_preset = _safe_part(preset) if preset else "default"
    tgt_file = _safe_part(Path(target_fasta).stem)
    qry_file = _safe_part(Path(query_fasta).stem)
    dirname = f"{qry_file}.{safe_query}_vs_{tgt_file}.{safe_target}"
    if output_dir is None:
        project_root = Path(__file__).resolve().parent.parent
        base_dir = project_root / "resources"
    else:
        base_dir = Path(output_dir)
    folder = base_dir / dirname
    filename = f"{dirname}.{safe_preset}.paf"
    return folder / filename


def _materialize_single_sequence_fasta(
    fasta_path: Path, sequence_name: str, work_dir: Path
) -> Path:
    """Write a tiny FASTA containing only one sequence for faster alignment."""
    sequences = read_fasta_sequences(fasta_path, select_names={sequence_name})
    if sequence_name not in sequences:
        raise ValueError(f"Sequence '{sequence_name}' not found in {fasta_path}")
    out_path = work_dir / f"{sequence_name}.fa"
    write_fasta({sequence_name: sequences[sequence_name]}, out_path)
    return out_path


def run_minimap2_alignment(
    target_fasta: Path,
    query_fasta: Path,
    target_seq: str,
    query_seq: str,
    *,
    output_path: Optional[Path] = None,
    preset: str = "asm10",
    threads: int = 4,
    minimap2: str = "minimap2",
    extra_args: Optional[Sequence[str]] = None,
    reuse_existing: bool = True,
    filter_sequences: bool = True,
    dry_run: bool = False,
    reverse_query: bool = False,
) -> AlignmentRun:
    """
    Run minimap2 (via mappy) to align `query_seq` onto `target_seq` and save a PAF file.

    Parameters mirror minimap2 wherever possible:
    - preset: passed to `-x`
    - threads: passed to `-t`
    - extra_args: kept for API compatibility; currently not forwarded in mappy mode

    File handling niceties:
    - If output_path is None or a directory, a folder named `{query}_vs_{target}`
      is created in CWD and a PAF of the same name is written inside.
    - When reuse_existing is True and the PAF already exists, the command is skipped.
    - Only the requested sequences are materialized (filter_sequences=True) to avoid
      full-assembly runtime. Disable if you want the original FASTA untouched.
    - reverse_query=True will reverse-complement the query before alignment.
    """
    output_path = (
        default_paf_path(
            target_fasta, query_fasta, target_seq, query_seq, preset, output_path, reverse_query=reverse_query
        )
        if output_path is None or Path(output_path).is_dir()
        else Path(output_path)
    )
    output_path = output_path.with_suffix(".paf")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if reuse_existing and output_path.exists():
        return AlignmentRun(
            target_fasta=Path(target_fasta),
            query_fasta=Path(query_fasta),
            target_seq=target_seq,
            query_seq=query_seq,
            preset=preset,
            output_path=output_path,
            cmd=[],
            skipped=True,
        )

    work_dir = Path(tempfile.mkdtemp(prefix="gapneedle_", dir=output_path.parent))
    try:
        target_records = read_fasta_sequences(target_fasta, select_names={target_seq})
        if target_seq not in target_records:
            raise ValueError(f"Sequence '{target_seq}' not found in {target_fasta}")
        t_seq_str = target_records[target_seq]
        t_len = len(t_seq_str)

        if filter_sequences:
            tgt_fa = work_dir / f"{target_seq}.fa"
            write_fasta({target_seq: t_seq_str}, tgt_fa)
        else:
            tgt_fa = Path(target_fasta)

        query_records = read_fasta_sequences(query_fasta, select_names={query_seq})
        if query_seq not in query_records:
            raise ValueError(f"Sequence '{query_seq}' not found in {query_fasta}")
        qry_seq = query_records[query_seq]
        if reverse_query:
            qry_seq = reverse_complement(qry_seq)
        q_len = len(qry_seq)

        cmd: List[str] = [f"mappy.Aligner({tgt_fa}, preset={preset}, n_threads={threads})"]
        run_info = AlignmentRun(
            target_fasta=Path(target_fasta),
            query_fasta=Path(query_fasta),
            target_seq=target_seq,
            query_seq=query_seq,
            preset=preset,
            output_path=output_path,
            cmd=cmd,
            skipped=False,
        )

        if dry_run:
            return run_info

        # mappy.Aligner expects reference file as positional arg
        aligner = mp.Aligner(str(tgt_fa), preset=preset, n_threads=threads)
        if not aligner:
            raise Minimap2Error("无法初始化 mappy.Aligner，请检查参考序列或参数。")

        with output_path.open("w") as out:
            # mappy map API does not accept seq_name; we propagate query_seq manually.
            for hit in aligner.map(qry_seq):
                # ctg may be missing/NA when aligning in-memory; fall back to target_seq
                t_name = getattr(hit, "ctg", None) or target_seq
                strand = "+" if hit.strand >= 0 else "-"
                t_len_hit = getattr(hit, "ctg_len", t_len)
                paf_fields = [
                    query_seq,          # qname
                    q_len,              # qlen
                    hit.q_st,
                    hit.q_en,
                    strand,
                    t_name,             # tname
                    t_len_hit,          # tlen
                    hit.r_st,
                    hit.r_en,
                    getattr(hit, "mlen", hit.blen),
                    hit.blen,
                    hit.mapq,
                ]
                out.write("\t".join(map(str, paf_fields)))
                cigar = getattr(hit, "cigar_str", None)
                if cigar:
                    out.write(f"\tcg:Z:{cigar}")
                if getattr(hit, "cs", None):
                    out.write(f"\tcs:Z:{hit.cs}")
                out.write("\n")

        return run_info
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
