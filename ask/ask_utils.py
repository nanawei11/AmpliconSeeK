from __future__ import annotations

import time
from pathlib import Path
from typing import Any


MODULE_DIR = Path(__file__).resolve().parent
DATA_DIR = MODULE_DIR.parent / "data"


def ask_output_paths(outprefix: str | Path) -> dict[str, Path]:
    prefix = Path(outprefix)
    return {
        "outfig": Path(str(prefix) + "_ask_plot"),
        "output_stats": Path(str(prefix) + "_ask_stats.tsv"),
        "output_cnamp": Path(str(prefix) + "_ask_amplified_segment.tsv"),
        "output_cnseg": Path(str(prefix) + "_ask_cn_segmentation.tsv"),
        "output_bpcand": Path(str(prefix) + "_ask_breakpoint.tsv"),
        "output_bpduo": Path(str(prefix) + "_ask_breakpoint_pair_raw.tsv"),
        "output_bppair": Path(str(prefix) + "_ask_breakpoint_pair.tsv"),
        "output_seg": Path(str(prefix) + "_ask_breakpoint_seg.tsv"),
        "output_circ": Path(str(prefix) + "_ask_amplicon_circular.tsv"),
        "output_line": Path(str(prefix) + "_ask_amplicon_linear.tsv"),
        "output_clip": Path(str(prefix) + "_ask_clip_count.bedgraph"),
        "output_bincnt": Path(str(prefix) + "_ask_bin_count.tsv"),
        "output_binnorm": Path(str(prefix) + "_ask_bin_count_norm.tsv"),
        "output_align": Path(str(prefix) + "_ask_alignment_sequence.tsv"),
        "output_align_dir": Path(str(prefix) + "_ask_junctionseq"),
        "output_circ_stat": Path(str(prefix) + "_ask_amplicon_circular_stat.tsv"),
        "output_jcs": Path(str(prefix) + "_ask_jcs.tsv"),
        "output_pdat_1": Path(str(prefix) + "_ask_step1.pdat"),
        "output_pdat_2": Path(str(prefix) + "_ask_step2.pdat"),
        "output_pdat_3": Path(str(prefix) + "_ask_step3.pdat"),
        "output_pdat_4": Path(str(prefix) + "_ask_step4.pdat"),
    }


def prepare_ask_output_paths(outprefix: str | Path) -> dict[str, Path]:
    paths = ask_output_paths(outprefix)
    Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
    paths["outfig"].mkdir(parents=True, exist_ok=True)
    paths["output_align_dir"].mkdir(parents=True, exist_ok=True)
    return paths


def first_existing_path(paths: list[Path | str | None]) -> Path | None:
    for path in paths:
        if path and Path(path).exists():
            return Path(path)
    return None


def resolve_ask_annotation_files(
    genome: str = "hg38",
    genefile: str | Path | None = None,
    sefile: str | Path | None = None,
    cgfile: str | Path | None = None,
    blfile: str | Path | None = None,
    gsfile: str | Path | None = None,
    biasfile: str | Path | None = None,
    data_dir: str | Path | None = None,
) -> dict[str, Path | None]:
    root = Path(data_dir) if data_dir else DATA_DIR
    return {
        "blfile": first_existing_path([blfile, root / f"{genome}_blacklist.bed"]),
        "genefile": first_existing_path(
            [
                genefile,
                root / f"{genome}_refgene_process.bed12",
                root / f"{genome}_refgene.bed12",
                root / f"{genome}_refgene.bed12.gz",
            ]
        ),
        "sefile": first_existing_path([sefile, root / f"se_{genome}_sort.bed"]),
        "gsfile": first_existing_path([gsfile, root / f"{genome}.genome"]),
        "cgfile": first_existing_path([cgfile, root / "Census_all_20200624_14_22_39.tsv"]),
        "biasfile": first_existing_path([biasfile, root / f"{genome}_bias.bed.gz"]),
    }


def plot_ask_amplicons(
    ask_module: Any,
    circ_anno: Any,
    line_anno: Any,
    cn_amp: Any,
    genefile: str | Path | None,
    sefile: str | Path | None,
    bin_norm: Any,
    fig_dir: str | Path,
    binsize: int = 10000,
    plot_n: int = 5,
    fig_width: int = 15,
    fontsize: int = 12,
    ext: float = 0.3,
    start_time: float | None = None,
    clean_old_pdf: bool = True,
) -> None:
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    if clean_old_pdf:
        for old_pdf in fig_dir.glob("*.pdf"):
            old_pdf.unlink()
    start = time.time() if start_time is None else start_time
    try:
        ask_module.plot_amplicon(
            circ_anno,
            line_anno,
            cn_amp,
            str(genefile) if genefile is not None else None,
            str(sefile) if sefile is not None else None,
            bin_norm,
            binsize=binsize,
            fig_dir=str(fig_dir),
            plot_n=plot_n,
            fig_width=fig_width,
            fontsize=fontsize,
            ext=ext,
        )
        print("plot amplicons - done             - {0:0.1f} seconds\n".format(time.time() - start))
    except Exception:
        print("plot amplicons - error             - {0:0.1f} seconds\n".format(time.time() - start))


def write_ask_stats(
    path: str | Path,
    cn_amp: Any,
    bp_cand_stats: Any,
    bp_duo: Any,
    circ_anno: Any,
    line_anno: Any,
) -> None:
    circ_n = circ_anno["AmpliconID"].nunique() if "AmpliconID" in circ_anno.columns else 0
    line_n = line_anno["AmpliconID"].nunique() if "AmpliconID" in line_anno.columns else 0
    with Path(path).open("w") as handle:
        handle.write(f"# of amplified segments: {len(cn_amp)}\n")
        handle.write(f"# of candidate breakpoints: {len(bp_cand_stats)}\n")
        handle.write(f"# of breakpoint pairs: {len(bp_duo)}\n")
        handle.write(f"# of circular amplicons: {circ_n}\n")
        handle.write(f"# of linear amplicons: {line_n}\n")
