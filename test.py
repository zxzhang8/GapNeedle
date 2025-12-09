from pathlib import Path

from gapfillet import GapFillet, manual_stitch_by_coordinates
from gapfillet.stitcher import telomere_details

if __name__ == '__main__':

    gf = GapFillet()
    tgt = Path("/home/zxzhang/remote/mouse/ICR/ICR.ont.p_ctg.fa")
    qry = Path("/home/zxzhang/remote/mouse/ICR/ICR.all.p_ctg.fa")
    tgt_seq = "ptg000015l"
    qry_seq = "ptg000001l"

    # run = gf.align(
    #     target_fasta=Path("/home/zxzhang/remote/mouse/ICR/ICR.ont.p_ctg.fa"),
    #     query_fasta=Path("/home/zxzhang/remote/mouse/ICR/ICR.all.p_ctg.fa"),
    #     target_seq="ptg000015l",
    #     query_seq="ptg000001l",
    #     threads=16,
    #     preset="asm5"
    # )
    #
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
    )

