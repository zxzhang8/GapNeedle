from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

from .aligner import (
    AlignmentRun,
    Minimap2Error,
    default_paf_path,
    run_minimap2_alignment,
)
from .io import read_fasta_sequences, reverse_complement, write_fasta


MessageFn = Callable[[str], None]
TELOMERE_MOTIF = "CCCTAA"
TELOMERE_WINDOW = 1000
TELOMERE_MIN_REPEATS = 15
HIGHLIGHT = "\x1b[1;33m"
OK_COLOR = "\x1b[1;32m"
RESET = "\x1b[0m"
MD_OK = '<span style="color:green;font-weight:bold">'
MD_WARN = '<span style="color:orange;font-weight:bold">'
MD_RESET = "</span>"


@dataclass
class PafRecord:
    target: str
    t_len: int
    t_start: int
    t_end: int
    strand: str
    query: str
    q_len: int
    q_start: int
    q_end: int
    matches: int
    aln_len: int
    mapq: int
    extras: List[str]

    @property
    def overlap(self) -> int:
        return min(self.t_end - self.t_start, self.q_end - self.q_start)


def parse_paf(path: Path, target_seq: str, query_seq: str) -> List[PafRecord]:
    records: List[PafRecord] = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            q_name, q_len, q_start, q_end, strand, t_name, t_len, t_start, t_end, matches, aln_len, mapq, *rest = parts
            if t_name != target_seq or q_name != query_seq:
                continue
            records.append(
                PafRecord(
                    target=t_name,
                    t_len=int(t_len),
                    t_start=int(t_start),
                    t_end=int(t_end),
                    strand=strand,
                    query=q_name,
                    q_len=int(q_len),
                    q_start=int(q_start),
                    q_end=int(q_end),
                    matches=int(matches),
                    aln_len=int(aln_len),
                    mapq=int(mapq),
                    extras=rest,
                )
            )
    return records


def suggest_overlaps(
    paf_path: Path,
    target_seq: str,
    query_seq: str,
    limit: Optional[int] = 10,
) -> List[PafRecord]:
    """Return candidate overlaps ordered by longest first."""
    records = parse_paf(paf_path, target_seq, query_seq)
    records.sort(key=lambda r: r.overlap, reverse=True)
    if limit:
        return records[:limit]
    return records


def _show_candidates(records: Sequence[PafRecord], say: MessageFn) -> None:
    say("可用的拼接候选（按overlap长度排序）:")
    for idx, rec in enumerate(records):
        q_len = f"{rec.q_len:,}"
        t_len = f"{rec.t_len:,}"
        overlap = f"{rec.overlap:,}"
        t_span = f"{rec.t_start:,}-{rec.t_end:,}"
        q_span = f"{rec.q_start:,}-{rec.q_end:,}"
        say(
            f"[{idx}] {rec.query}({q_len}bp) -> {rec.target}({t_len}bp) "
            f"strand={rec.strand} overlap={overlap} mapq={rec.mapq} "
            f"t:{t_span} q:{q_span}"
        )


def _orientation_adjustment(rec: PafRecord, seq: str) -> Tuple[str, int, int]:
    """Return oriented query sequence and adjusted start/end on that orientation."""
    if rec.strand == "+":
        return seq, rec.q_start, rec.q_end
    oriented = reverse_complement(seq)
    return oriented, rec.q_len - rec.q_end, rec.q_len - rec.q_start


def _junction_preview(
    left: str, right: str, context: int, highlight_len: int = 10, color: str = HIGHLIGHT
) -> str:
    """
    Build a single-line preview with colored junction:
    ...<left tail><color>|<right head><reset>...
    """
    ctx = max(0, context)
    l_tail = left[-ctx:]
    r_head = right[:ctx]
    hl = min(highlight_len, len(l_tail), len(r_head))
    l_head_plain = l_tail[: len(l_tail) - hl]
    l_hl = l_tail[len(l_tail) - hl :]
    r_hl = r_head[:hl]
    r_rest = r_head[hl:]
    if hl == 0:
        return f"{l_tail}|{r_head}"
    return f"{l_head_plain}{color}{l_hl}|{r_hl}{RESET}{r_rest}"


def _highlight_diff(a: str, b: str, color: str = HIGHLIGHT, reset: str = RESET) -> Tuple[str, str]:
    """Return two strings with differing positions wrapped in color."""
    res_a: List[str] = []
    res_b: List[str] = []
    for x, y in zip(a, b):
        if x == y:
            res_a.append(x)
            res_b.append(y)
        else:
            res_a.append(f"{color}{x}{reset}")
            res_b.append(f"{color}{y}{reset}")
    return "".join(res_a), "".join(res_b)


def _md_wrap(text: str, color: str) -> str:
    return f'<span style="color:{color};font-weight:bold">{text}</span>'


def _highlight_diff_md(a: str, b: str) -> Tuple[str, str]:
    """Markdown-friendly diff highlight (orange)."""
    return _highlight_diff(a, b, color=MD_WARN, reset=MD_RESET)


def _junction_preview_md(left: str, right: str, context: int, highlight_len: int = 10) -> str:
    """
    Markdown preview with colored junction.
    区域一致则绿色，否 orange。
    """
    ctx = max(0, context)
    l_tail = left[-ctx:]
    r_head = right[:ctx]
    hl = min(highlight_len, len(l_tail), len(r_head))
    if hl == 0:
        return f"{l_tail}|{r_head}"
    l_head_plain = l_tail[: len(l_tail) - hl]
    l_hl = l_tail[len(l_tail) - hl :]
    r_hl = r_head[:hl]
    r_rest = r_head[hl:]
    color = "green" if l_hl == r_hl else "orange"
    return f"{l_head_plain}{_md_wrap(l_hl, color)}|{_md_wrap(r_hl, color)}{r_rest}"


def _find_runs(subseq: str, motif: str) -> List[Tuple[int, int, List[Tuple[int, int]]]]:
    """Find all consecutive runs of motif in subseq; returns (count, span_end, matches)."""
    seq = subseq.upper()
    motif = motif.upper()
    mlen = len(motif)
    runs: List[Tuple[int, int, List[Tuple[int, int]]]] = []
    i = 0
    while True:
        pos = seq.find(motif, i)
        if pos == -1:
            break
        matches: List[Tuple[int, int]] = []
        j = pos
        while seq.startswith(motif, j):
            matches.append((j, j + mlen))
            j += mlen
        runs.append((len(matches), j, matches))
        i = pos + 1
    return runs


def _continuous_telomere_details(subseq: str, motif: str, min_repeats: int):
    """
    Check for consecutive telomere motif repeats (or reverse complement).
    Returns (has, matches, span_len, mutated_bases, repeat_count).
    Only runs with repeat_count >= min_repeats are considered hits.
    """
    motif = motif.upper()
    motif_rc = reverse_complement(motif)
    runs = _find_runs(subseq, motif) + _find_runs(subseq, motif_rc)
    if not runs:
        return False, [], 0, 0, 0
    # pick the longest run
    best = max(runs, key=lambda x: x[0])
    count, span_end, matches = best
    if count < min_repeats:
        return False, [], 0, 0, 0
    span_start = matches[0][0]
    span_len = span_end - span_start
    mutated = max(span_len - count * len(motif), 0)  # should be 0 for perfect repeats
    return True, matches, span_len, mutated, count


def _telomere_flags(
    seq: str, window: int = TELOMERE_WINDOW, min_repeats: int = TELOMERE_MIN_REPEATS
) -> Tuple[bool, bool]:
    """Lightweight flag check on a raw sequence string."""
    window = min(window, len(seq))
    left_sub = seq[:window]
    right_sub = seq[-window:]
    left_has, *_ = _continuous_telomere_details(left_sub, TELOMERE_MOTIF, min_repeats)
    right_has, *_ = _continuous_telomere_details(right_sub, TELOMERE_MOTIF, min_repeats)
    return left_has, right_has


def telomere_details(
    fasta_path: Path,
    seq_name: str,
    *,
    window: int = 10_000_000,
    motif: str = "CCCTAA",
    min_repeats: int = TELOMERE_MIN_REPEATS,
) -> Tuple[dict, dict]:
    """
    Return telomere motif details (consecutive repeats) for both ends.
    Each side dict: {has, matches, span_len, mutated_bases, repeat_count}
    - matches: list of (start, end) within the window (0-based) for the best run
    """
    seqs = read_fasta_sequences(fasta_path, select_names={seq_name})
    if seq_name not in seqs:
        raise ValueError(f"Sequence '{seq_name}' not found in {fasta_path}")
    seq = seqs[seq_name]
    window = min(window, len(seq))
    motif = motif.upper()

    left_sub = seq[:window]
    right_sub = seq[-window:]
    left_has, left_matches, left_span, left_mut, left_count = _continuous_telomere_details(
        left_sub, motif, min_repeats
    )
    right_has, right_matches, right_span, right_mut, right_count = _continuous_telomere_details(
        right_sub, motif, min_repeats
    )
    left_info = {
        "has": left_has,
        "matches": left_matches,
        "span_len": left_span,
        "mutated_bases": left_mut,
        "repeat_count": left_count,
    }
    right_info = {
        "has": right_has,
        "matches": right_matches,
        "span_len": right_span,
        "mutated_bases": right_mut,
        "repeat_count": right_count,
    }
    return left_info, right_info


def telomere_presence(
    fasta_path: Path,
    seq_name: str,
    *,
    window: int = 10_000_000,
    motif: str = "CCCTAA",
    min_repeats: int = TELOMERE_MIN_REPEATS,
) -> Tuple[bool, bool]:
    """
    Check whether a sequence contains consecutive telomere motif repeats
    (or its reverse complement) within `window` bases from both ends.

    Returns (has_left, has_right).
    """
    left_info, right_info = telomere_details(
        fasta_path, seq_name, window=window, motif=motif, min_repeats=min_repeats
    )
    return left_info["has"], right_info["has"]


def manual_stitch_by_coordinates(
    target_fasta: Path,
    query_fasta: Path,
    target_seq: str,
    query_seq: str,
    *,
    say: MessageFn = print,
    context: int = 200,
    stop_token: str = "x",
    reverse_query: bool = False,
) -> str:
    """
    Interactively stitch segments from two sequences by coordinates.

    - User chooses which sequence (target/query) and provides start-end (0-based, end-exclusive).
    - Prints ~context bp around both breakpoints: left breakpoint (before/after start),
      then right breakpoint (before/after end). 对于第 2 段起，先展示左断点，再在加入后展示右断点。
    - After each addition, compares previous segment tail vs current head to hint consistency.
    - Stops when user inputs stop_token (default 'x') as the source.
    - At the end, prints all junction contexts with a colored '|' marker, then asks for FASTA path/name to save.
    - reverse_query=True 将查询序列反向互补后再用于拼接（坐标按反向后的序列计算）。
    - Returns the merged sequence string (regardless of saving).
    """
    t_records = read_fasta_sequences(target_fasta, select_names={target_seq})
    q_records = read_fasta_sequences(query_fasta, select_names={query_seq})
    if target_seq not in t_records or query_seq not in q_records:
        raise ValueError("FASTA中找不到目标序列或查询序列。")
    q_seq = q_records[query_seq]
    if reverse_query:
        q_seq = reverse_complement(q_seq)
        say("已使用查询序列的反向互补进行手动拼接。")
    seq_map = {"t": (target_seq, t_records[target_seq]), "q": (query_seq, q_seq)}
    stop_token = stop_token.lower()
    if stop_token in seq_map:
        raise ValueError(f"stop_token '{stop_token}' 不能与来源键重复，请使用其他字符。")

    say(
        f"手动拼接模式：输入来源(t/q) 与坐标 start-end（0-based，end 不含），"
        f"输入 {stop_token} 退出。"
    )
    pieces: List[Tuple[str, int, int, str, str, str, str, str]] = []
    merged_parts: List[str] = []
    while True:
        src = input(f"选择来源 [t={target_seq}/q={query_seq}/{stop_token}=退出]: ").strip().lower()
        if src == stop_token:
            break
        if src not in seq_map:
            say("无效来源，请输入 t/q 或退出指令。")
            continue
        coord_raw = input("输入坐标 start-end (0-based, 末端不含): ").strip()
        if "-" not in coord_raw:
            say("坐标格式应为 start-end。")
            continue
        try:
            start_s, end_s = coord_raw.split("-", 1)
            start = int(start_s)
            end = int(end_s)
        except ValueError:
            say("坐标必须是整数。")
            continue
        name, seq = seq_map[src]
        if not (0 <= start < end <= len(seq)):
            say(f"坐标越界，序列长度 {len(seq):,}。")
            continue
        seg = seq[start:end]
        left_before = seq[max(0, start - context) : start]
        left_after = seq[start : min(len(seq), start + context)]
        right_before = seq[max(0, end - context) : end]
        right_after = seq[end : min(len(seq), end + context)]

        say(f"{name}:{start}-{end} 长度{len(seg):,}bp")
        say(
            f"左断点 前{len(left_before)}bp: {left_before}\n"
            f"左断点 后{len(left_after)}bp: {left_after}"
        )

        if pieces:
            # Compare breakpoint contexts: left vs left, right vs right.
            prev_left = pieces[-1][6]  # 上一段断点左侧（末端前 context）
            prev_right = pieces[-1][7]  # 上一段断点右侧（末端后 context）
            curr_left = left_before  # 当前段断点左侧
            curr_right = left_after  # 当前段断点右侧

            say("与上一段断点一致性检查:")

            left_len = min(50, len(prev_left), len(curr_left))
            if left_len:
                prev_left_slice = prev_left[-left_len:]
                curr_left_slice = curr_left[-left_len:]
                if prev_left_slice == curr_left_slice:
                    say(f"左侧一致（末{left_len}bp）: {OK_COLOR}{prev_left_slice}{RESET}")
                else:
                    hi_prev, hi_curr = _highlight_diff(prev_left_slice, curr_left_slice)
                    say(f"左侧不一致（末{left_len}bp，高亮为差异位点）:")
                    say(f"上一段左: {hi_prev}")
                    say(f"当前段左: {hi_curr}")
            else:
                say("左侧无可比对的上下文。")

            right_len = min(50, len(prev_right), len(curr_right))
            if right_len:
                prev_right_slice = prev_right[:right_len]
                curr_right_slice = curr_right[:right_len]
                if prev_right_slice == curr_right_slice:
                    say(f"右侧一致（首{right_len}bp）: {OK_COLOR}{prev_right_slice}{RESET}")
                else:
                    hi_prev, hi_curr = _highlight_diff(prev_right_slice, curr_right_slice)
                    say(f"右侧不一致（首{right_len}bp，高亮为差异位点）:")
                    say(f"上一段右: {hi_prev}")
                    say(f"当前段右: {hi_curr}")
            else:
                say("右侧无可比对的上下文。")

        pieces.append(
            (
                src,
                start,
                end,
                seg,
                left_before,
                left_after,
                right_before,
                right_after,
            )
        )
        merged_parts.append(seg)

        say(
            f"右断点 前{len(right_before)}bp: {right_before}\n"
            f"右断点 后{len(right_after)}bp: {right_after}"
        )

    merged_seq = "".join(merged_parts)
    say(f"手动拼接完成，共 {len(pieces)} 段，合并长度 {len(merged_seq):,}bp")
    if len(pieces) >= 2:
        say(f"断点汇总（上下文 {context}bp，竖线为断点，左右高亮）:")
        for i in range(len(pieces) - 1):
            l_src, l_start, l_end, l_seg, _, _, l_right_before, _ = pieces[i]
            r_src, r_start, r_end, r_seg, r_left_before, r_left_after, _, _ = pieces[i + 1]
            preview = _junction_preview(l_right_before, r_left_after, context)
            say(
                f"[{i}] {l_src}:{l_start}-{l_end} -> {r_src}:{r_start}-{r_end}\n"
                f"{preview}"
            )

    save_ans = input("是否保存为 FASTA? [y/N]: ").strip().lower()
    if save_ans in {"y", "yes"}:
        while True:
            out_path_raw = input("输出文件路径: ").strip()
            if out_path_raw:
                out_path = Path(out_path_raw)
                break
            say("路径不能为空。")
        default_name = f"{target_seq}+{query_seq}"
        seq_name = input(f"输出序列名 (默认 {default_name}): ").strip() or default_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta({seq_name: merged_seq}, out_path)
        log_path = out_path.with_suffix(".md")
        log_lines = [
            "# 手动拼接记录",
            "",
            f"- 目标 FASTA: {target_fasta}",
            f"- 查询 FASTA: {query_fasta}",
            f"- 目标序列: {target_seq}",
            f"- 查询序列: {query_seq}",
            f"- 输出 FASTA: {out_path}",
            f"- 输出序列名: {seq_name}",
            f"- 段数: {len(pieces)}",
            f"- 合并长度: {len(merged_seq):,} bp",
            "",
            "## 选取的片段",
        ]
        for idx, (src, start, end, seg, _, _, _, _) in enumerate(pieces):
            src_label = "目标" if src == "t" else "查询"
            log_lines.append(f"- [{idx}] {src_label}({src}) {start}-{end} 长度 {len(seg):,} bp")
        if len(pieces) >= 2:
            log_lines.append("")
            log_lines.append(f"## 断点上下文（各取 {context}bp）")
            for i in range(len(pieces) - 1):
                (
                    l_src,
                    l_start,
                    l_end,
                    _,
                    _,
                    _,
                    l_right_before,
                    l_right_after,
                ) = pieces[i]
                (
                    r_src,
                    r_start,
                    r_end,
                    _,
                    r_left_before,
                    r_left_after,
                    _,
                    _,
                ) = pieces[i + 1]

                preview = _junction_preview_md(l_right_before, r_left_after, context)
                log_lines.append(f"- [{i}] {l_src}:{l_start}-{l_end} -> {r_src}:{r_start}-{r_end}")
                log_lines.append(f"  - 断点上下文: {preview}")

                left_len = min(50, len(l_right_before), len(r_left_before))
                if left_len:
                    prev_slice = l_right_before[-left_len:]
                    curr_slice = r_left_before[-left_len:]
                    if prev_slice == curr_slice:
                        log_lines.append(
                            f"  - 左侧一致（末{left_len}bp）: {_md_wrap(prev_slice, 'green')}"
                        )
                    else:
                        hi_prev, hi_curr = _highlight_diff_md(prev_slice, curr_slice)
                        log_lines.append(f"  - 左侧不一致（末{left_len}bp，高亮为差异位点）:")
                        log_lines.append(f"    上一段左: {hi_prev}")
                        log_lines.append(f"    当前段左: {hi_curr}")
                else:
                    log_lines.append("  - 左侧无可比对的上下文。")

                right_len = min(50, len(l_right_after), len(r_left_after))
                if right_len:
                    prev_slice = l_right_after[:right_len]
                    curr_slice = r_left_after[:right_len]
                    if prev_slice == curr_slice:
                        log_lines.append(
                            f"  - 右侧一致（首{right_len}bp）: {_md_wrap(prev_slice, 'green')}"
                        )
                    else:
                        hi_prev, hi_curr = _highlight_diff_md(prev_slice, curr_slice)
                        log_lines.append(f"  - 右侧不一致（首{right_len}bp，高亮为差异位点）:")
                        log_lines.append(f"    上一段右: {hi_prev}")
                        log_lines.append(f"    当前段右: {hi_curr}")
                else:
                    log_lines.append("  - 右侧无可比对的上下文。")
        log_path.write_text("\n".join(log_lines))
        say(f"拼接结果已保存到 {out_path} ，序列名 {seq_name}")
        say(f"操作记录已保存到 {log_path}")

    return merged_seq


def _choose_segment_source(
    name: str,
    options: List[Tuple[str, str, str]],
    default_key: str,
    interactive: bool,
    say: MessageFn,
) -> Tuple[str, str]:
    """
    Let users choose which sequence to use for a segment.

    options: list of (key, label, seq)
    returns: (chosen_key, chosen_seq)
    """
    available = [(k, label, seq) for k, label, seq in options if seq is not None]
    if not available:
        return default_key, ""

    keys = [k for k, _, _ in available]
    default = default_key if default_key in keys else keys[0]
    if not interactive:
        chosen_key, chosen_label, chosen_seq = next(
            (k, label, seq) for k, label, seq in available if k == default
        )
        return chosen_key, chosen_seq

    say(
        f"{name} 段可选来源: "
        + ", ".join(f"{k}={label}({len(seq):,}bp)" for k, label, seq in available)
    )
    while True:
        raw = input(
            f"选择 {name} 段来源 [{'/'.join(keys)}] (默认 {default}): "
        ).strip()
        if not raw:
            choice = default
        else:
            choice = raw.lower()
        if choice in keys:
            chosen_key, chosen_label, chosen_seq = next(
                (k, label, seq) for k, label, seq in available if k == choice
            )
            say(f"{name} 段使用 {chosen_label} 序列。")
            return chosen_key, chosen_seq
        say("输入无效，请输入可选的键。")


def stitch_from_paf(
    paf_path: Path,
    target_fasta: Path,
    query_fasta: Path,
    target_seq: str,
    query_seq: str,
    *,
    selection: Optional[int] = None,
    interactive: bool = True,
    confirm: bool = True,
    say: MessageFn = print,
    output_fasta: Optional[Path] = None,
) -> Path:
    """
    Stitch two sequences using an alignment chosen from a PAF.

    - If selection is None and interactive=True, users will be prompted to choose.
    - Resulting FASTA contains a single merged sequence named `{target}+{query}`.
    """
    paf_path = Path(paf_path)
    if not paf_path.exists():
        raise FileNotFoundError(f"PAF not found: {paf_path}")

    candidates = suggest_overlaps(paf_path, target_seq, query_seq)
    if not candidates:
        raise ValueError("PAF中未找到对应序列的比对结果，无法拼接。")

    target_seq_dict = read_fasta_sequences(target_fasta, select_names={target_seq})
    query_seq_dict = read_fasta_sequences(query_fasta, select_names={query_seq})
    if target_seq not in target_seq_dict or query_seq not in query_seq_dict:
        raise ValueError("FASTA中找不到目标序列或查询序列。")

    t_seq = target_seq_dict[target_seq]
    q_seq_raw = query_seq_dict[query_seq]
    t_tel_left, t_tel_right = _telomere_flags(t_seq)
    q_tel_left, q_tel_right = _telomere_flags(q_seq_raw)
    flag = lambda x: "是" if x else "否"
    say(
        f"端粒检测 (窗口{TELOMERE_WINDOW}bp，基序 {TELOMERE_MOTIF}/互补，最少连续 {TELOMERE_MIN_REPEATS} 次): "
        f"\n{target_seq} 左侧[{flag(t_tel_left)}] 右侧[{flag(t_tel_right)}]"
        f"\n{query_seq} 左侧[{flag(q_tel_left)}] 右侧[{flag(q_tel_right)}]"
    )

    if selection is None and interactive:
        _show_candidates(candidates, say)
        while True:
            raw = input("选择需要的拼接候选编号: ").strip()
            if raw.isdigit() and int(raw) < len(candidates):
                selection = int(raw)
                break
            say("输入无效，请输入列表中的编号。")
    elif selection is None:
        selection = 0

    rec = candidates[selection]
    identity = rec.matches / rec.aln_len if rec.aln_len else 0
    overlap = f"{rec.overlap:,}"
    t_span = f"{rec.t_start:,}-{rec.t_end:,}"
    q_span = f"{rec.q_start:,}-{rec.q_end:,}"
    say(
        f"已选择 [{selection}] strand={rec.strand} overlap={overlap} "
        f"t:{t_span} q:{q_span} "
        f"identity≈{identity:.3f}"
    )
    if interactive and confirm:
        ans = input("是否继续拼接? [y/N]: ").strip().lower()
        if ans not in {"y", "yes"}:
            raise KeyboardInterrupt("用户取消拼接。")

    q_seq, q_start, q_end = _orientation_adjustment(rec, q_seq_raw)

    # 计算可供选择的三段：左侧、重叠、右侧；额外的 overhang 自动拼上。
    left_common_len = min(rec.t_start, q_start)
    right_common_len = min(len(t_seq) - rec.t_end, len(q_seq) - q_end)
    overlap_len = rec.overlap

    left_common_opts = [
        ("t", "目标", t_seq[:left_common_len]),
        ("q", "查询", q_seq[:left_common_len]),
    ]
    overlap_opts = [
        ("t", "目标", t_seq[rec.t_start : rec.t_start + overlap_len]),
        ("q", "查询", q_seq[q_start : q_start + overlap_len]),
    ]
    right_common_opts = [
        ("t", "目标", t_seq[rec.t_end : rec.t_end + right_common_len]),
        ("q", "查询", q_seq[q_end : q_end + right_common_len]),
    ]

    left_overhang = (
        t_seq[left_common_len : rec.t_start]
        if rec.t_start > q_start
        else q_seq[left_common_len:q_start]
        if q_start > rec.t_start
        else ""
    )
    right_overhang = (
        t_seq[rec.t_end + right_common_len :]
        if len(t_seq) - rec.t_end > len(q_seq) - q_end
        else q_seq[q_end + right_common_len :]
        if len(q_seq) - q_end > len(t_seq) - rec.t_end
        else ""
    )

    left_choice, left_seq = _choose_segment_source(
        "左侧", left_common_opts, default_key="t", interactive=interactive, say=say
    )
    overlap_choice, overlap_seq = _choose_segment_source(
        "重叠区", overlap_opts, default_key="t", interactive=interactive, say=say
    )
    right_choice, right_seq = _choose_segment_source(
        "右侧", right_common_opts, default_key="t", interactive=interactive, say=say
    )

    merged_seq = "".join(
        [left_seq, left_overhang, overlap_seq, right_seq, right_overhang]
    )
    say(
        f"已选择 左侧={left_choice}({len(left_seq):,}bp) "
        f"重叠区={overlap_choice}({len(overlap_seq):,}bp) "
        f"右侧={right_choice}({len(right_seq):,}bp)，"
        f"自动附加 overhang 左 {len(left_overhang):,}bp 右 {len(right_overhang):,}bp，"
        f"合并后长度 {len(merged_seq):,}bp"
    )

    out_path = (
        Path(output_fasta)
        if output_fasta
        else paf_path.with_suffix(".stitched.fasta")
    )
    merged_name = f"{rec.target}+{rec.query}"
    write_fasta({merged_name: merged_seq}, out_path)
    say(f"拼接完成，结果已写入: {out_path}")
    return out_path


class GapFillet:
    """
    High-level helper that wraps alignment and stitching for Jupyter-friendly use.
    """

    def __init__(self, say: MessageFn = print):
        self.say = say

    def _print_paf_records(self, paf_path: Path) -> None:
        """Print all alignment records from a PAF file."""
        if not Path(paf_path).exists():
            self.say(f"PAF 文件不存在，无法打印比对记录: {paf_path}")
            return
        lines = Path(paf_path).read_text().splitlines()
        if not lines:
            self.say("PAF 为空，暂无比对记录。")
            return
        self.say("比对记录：")
        for idx, line in enumerate(lines):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                self.say(f"[{idx}] {line}")
                continue
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, matches, aln_len, mapq, *_ = parts
            arrow = "->" if strand == "+" else "<-"
            desc = (
                f"[{idx}]\t{qname}({int(qlen):,}bp)\t{int(qstart):,} - {int(qend):,}\t{arrow}\t"
                f"{tname}({int(tlen):,}bp)\t{int(tstart):,} - {int(tend):,}\t"
                f"matches={int(matches):,}/{int(aln_len):,} mapq={mapq}"
            )
            self.say(desc)

    def align(
        self,
        target_fasta: Path,
        query_fasta: Path,
        target_seq: str,
        query_seq: str,
        reverse_query: bool = False,
        **kwargs,
    ) -> AlignmentRun:
        self.say("调用 minimap2 生成 PAF ...")
        run = run_minimap2_alignment(
            target_fasta=target_fasta,
            query_fasta=query_fasta,
            target_seq=target_seq,
            query_seq=query_seq,
            reverse_query=reverse_query,
            **kwargs,
        )
        if run.skipped:
            self.say(f"检测到已存在的 PAF，跳过运行: {run.output_path}")
        else:
            self.say(f"PAF 已生成: {run.output_path}")
            if reverse_query:
                self.say("注意：本次比对使用了查询序列的反向互补。")
        self._print_paf_records(run.output_path)
        return run

    def stitch(
        self,
        target_fasta: Optional[Path] = None,
        query_fasta: Optional[Path] = None,
        target_seq: Optional[str] = None,
        query_seq: Optional[str] = None,
        *,
        paf_path: Optional[Path] = None,
        auto_align: bool = True,
        run: Optional[AlignmentRun] = None,
        preset: Optional[str] = None,
        reverse_query: bool = False,
        **kwargs,
    ) -> Path:
        """
        Stitch two sequences. If paf_path is missing and auto_align=True, minimap2 will run first.
        kwargs are forwarded to stitch_from_paf.
        """
        # Prefer explicit args; otherwise fall back to AlignmentRun.
        target_fasta = target_fasta or (run.target_fasta if run else None)
        query_fasta = query_fasta or (run.query_fasta if run else None)
        target_seq = target_seq or (run.target_seq if run else None)
        query_seq = query_seq or (run.query_seq if run else None)

        if not all([target_fasta, query_fasta, target_seq, query_seq]):
            raise ValueError("缺少 target/query fasta 或序列名；请传入参数或提供 run 对象。")

        preset = preset or getattr(run, "preset", "asm10")
        paf_path = paf_path or getattr(run, "output_path", None) or kwargs.pop("output_path", None)
        if paf_path is None:
            paf_path = default_paf_path(target_seq, query_seq, preset, reverse_query=reverse_query)

        if not Path(paf_path).exists():
            if not auto_align:
                raise FileNotFoundError(
                    f"PAF不存在 ({paf_path})，请先运行 align 或设置 auto_align=True。"
                )
            self.align(
                target_fasta=target_fasta,
                query_fasta=query_fasta,
                target_seq=target_seq,
                query_seq=query_seq,
                preset=preset,
                reverse_query=reverse_query,
            )

        return stitch_from_paf(
            paf_path=paf_path,
            target_fasta=target_fasta,
            query_fasta=query_fasta,
            target_seq=target_seq,
            query_seq=query_seq,
            say=self.say,
            **kwargs,
        )

    def scan_gaps(self, fasta_path: Path, min_gap: int = 10) -> List[Tuple[str, int, int]]:
        """
        Return gap intervals (N stretches) for quick inspection.
        Each tuple: (seq_name, start, end). start/end are 0-based half-open.
        """
        seqs = read_fasta_sequences(fasta_path)
        gaps: List[Tuple[str, int, int]] = []
        for name, seq in seqs.items():
            start = None
            for i, base in enumerate(seq):
                if base.upper() == "N":
                    if start is None:
                        start = i
                else:
                    if start is not None and i - start >= min_gap:
                        gaps.append((name, start, i))
                    start = None
            if start is not None and len(seq) - start >= min_gap:
                gaps.append((name, start, len(seq)))
        self.say(f"检测到 {len(gaps)} 个 gap 区间（长度≥{min_gap}）")
        return gaps

    def check_telomere(
        self,
        fasta_path: Path,
        seq_name: str,
        *,
        window: int = 10_000_000,
        motif: str = "CCCTAA",
        min_repeats: int = TELOMERE_MIN_REPEATS,
    ) -> Tuple[bool, bool]:
        """
        Query whether telomere motif (and its reverse complement) appears near both ends.
        Prints detail stats and returns (has_left, has_right).
        """
        left_info, right_info = telomere_details(
            fasta_path=fasta_path,
            seq_name=seq_name,
            window=window,
            motif=motif,
            min_repeats=min_repeats,
        )
        flag = lambda x: "是" if x else "否"
        self.say(
            f"{seq_name}: 左侧{flag(left_info['has'])} "
            f"(重复{left_info['repeat_count']}, 长度{left_info['span_len']:,}bp, 变异碱基{left_info['mutated_bases']})；"
            f"右侧{flag(right_info['has'])} "
            f"(重复{right_info['repeat_count']}, 长度{right_info['span_len']:,}bp, 变异碱基{right_info['mutated_bases']}) "
            f"(窗口 {window:,}bp, 基序 {motif}, 最低重复 {min_repeats})"
        )
        return left_info["has"], right_info["has"]
