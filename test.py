from pathlib import Path

from gapneedle import GapFillet, manual_stitch_by_coordinates
from gapneedle.stitcher import telomere_details

if __name__ == "__main__":
    gf = GapFillet()
    tgt = Path("/home/zxzhang/remote/mouse/ICR/ICR.ont.p_ctg.fa")
    qry = Path("/home/zxzhang/remote/mouse/ICR/ICR.all.p_ctg.fa")
    tgt_seq = "ptg000008l"
    qry_seq = "ptg000021l"
    reverse_query = False

    run = gf.align(
       target_fasta=tgt,
       query_fasta=qry,
       target_seq=tgt_seq,
       query_seq=qry_seq,
       threads=16,
       preset="asm20",
        reverse_query=reverse_query,
    )

    # merged = gf.stitch(run=run)

    # left, right = telomere_details(
    #     Path("/home/zxzhang/remote/mouse/ICR/ICR.all.p_ctg.fa"),
    #     "ptg000013l",
    #     window=8_000_000,
    #     min_repeats=5
    # )
    # pass

    manual_stitch_by_coordinates(
        tgt,
        qry,
        tgt_seq,
        qry_seq,
        reverse_query=reverse_query,
    )
