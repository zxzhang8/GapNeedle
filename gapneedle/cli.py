"""
Minimal CLI wrapper so users can try GapNeedle without writing code.
"""
import argparse
from pathlib import Path

from .stitcher import GapNeedle


def main(argv=None):
    parser = argparse.ArgumentParser(description="GapNeedle: minimap2 + stitching helper")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_align = sub.add_parser("align", help="运行 minimap2 生成 PAF")
    p_align.add_argument("target_fasta")
    p_align.add_argument("query_fasta")
    p_align.add_argument("target_seq")
    p_align.add_argument("query_seq")
    p_align.add_argument("--preset", default="asm10")
    p_align.add_argument("--threads", type=int, default=4)
    p_align.add_argument("--output", type=Path, default=None, help="输出文件或目录")
    p_align.add_argument("--minimap2", default="minimap2")
    p_align.add_argument("--no-filter", action="store_true", help="不抽取单条序列")

    p_stitch = sub.add_parser("stitch", help="根据 PAF 拼接")
    p_stitch.add_argument("paf")
    p_stitch.add_argument("target_fasta")
    p_stitch.add_argument("query_fasta")
    p_stitch.add_argument("target_seq")
    p_stitch.add_argument("query_seq")
    p_stitch.add_argument("--selection", type=int, default=None, help="候选编号，不填则交互选择")
    p_stitch.add_argument("--output", type=Path, default=None, help="拼接结果 FASTA 路径")

    p_scan = sub.add_parser("scan-gaps", help="扫描 N 区域")
    p_scan.add_argument("fasta")
    p_scan.add_argument("--min-gap", type=int, default=10)

    args = parser.parse_args(argv)
    gf = GapNeedle()

    if args.cmd == "align":
        run = gf.align(
            target_fasta=args.target_fasta,
            query_fasta=args.query_fasta,
            target_seq=args.target_seq,
            query_seq=args.query_seq,
            preset=args.preset,
            threads=args.threads,
            output_path=args.output,
            minimap2=args.minimap2,
            filter_sequences=not args.no_filter,
        )
        if not run.skipped:
            print("minimap2 命令:", " ".join(run.cmd))
        print("PAF:", run.output_path)
    elif args.cmd == "stitch":
        gf.stitch(
            paf_path=args.paf,
            target_fasta=args.target_fasta,
            query_fasta=args.query_fasta,
            target_seq=args.target_seq,
            query_seq=args.query_seq,
            selection=args.selection,
            interactive=args.selection is None,
            output_fasta=args.output,
        )
    elif args.cmd == "scan-gaps":
        gaps = gf.scan_gaps(args.fasta, min_gap=args.min_gap)
        for name, start, end in gaps:
            print(f"{name}\t{start}\t{end}")


if __name__ == "__main__":
    main()
