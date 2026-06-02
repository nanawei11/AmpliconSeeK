#!/usr/bin/env python3
"""
ASK-search

Search single-cell ATAC evidence for known ecDNA breakpoint junctions.

The module is designed around the workflow described in
ecDNA-Search Module Methodology:

1. Read ASK circular amplicon calls.
2. Convert each circular structure into ordered breakpoint pairs.
3. Match those breakpoint pairs to ASK breakpoint-pair evidence.
4. Map supporting read names to cell barcodes.
5. Build cell-by-breakpoint and cell-by-ecDNA evidence tables.
6. Classify ecDNA-positive cells by a configurable support threshold.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import itertools
import json
import pickle
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import sys

import pandas as pd
from pandas.errors import EmptyDataError
try:
    from .ask_utils import (
        DATA_DIR,
        ask_output_paths,
        plot_ask_amplicons,
        prepare_ask_output_paths,
        resolve_ask_annotation_files,
        write_ask_stats,
    )
except ImportError:  # Allow direct script execution: python ask/ecDNA_search.py ...
    from ask_utils import (
        DATA_DIR,
        ask_output_paths,
        plot_ask_amplicons,
        prepare_ask_output_paths,
        resolve_ask_annotation_files,
        write_ask_stats,
    )

try:
    import pysam
except ImportError:  # pragma: no cover - notebook users may install pysam later.
    pysam = None


TARGET_GENES = ("EGFR", "MDM4", "PDGFRA")

MODULE_DIR = Path(__file__).resolve().parent

DEFAULT_ASK_PATHS = (
    MODULE_DIR,
)

DEFAULT_ASK_DATA_PATHS = (
    DATA_DIR,
)


@dataclass(frozen=True)
class SampleInput:
    sample_id: str
    circular_path: Path | pd.DataFrame
    filename: str | None = None


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Input table does not exist: {path}. "
            "Please replace the placeholder path with your real ecDNA structure file."
        )
    if path.stat().st_size == 0:
        raise ValueError(
            f"Input table is empty: {path}. "
            "The known ecDNA structure file must contain columns such as "
            "AmpliconID, Chrom, Start, End, and Strand."
        )
    try:
        if path.suffix.lower() in {".csv"}:
            return pd.read_csv(path)
        return pd.read_csv(path, sep="\t")
    except EmptyDataError as exc:
        raise ValueError(
            f"No columns could be parsed from {path}. "
            "Check that the file is not empty and is a CSV/TSV table with a header."
        ) from exc


ASK_CIRCULAR_NOHEADER_COLUMNS = [
    "Chrom",
    "Start",
    "End",
    "Strand",
    "Score1",
    "Score2",
    "AmpliconID",
    "Gene",
    "CancerGene",
    "SegmentEdge",
]


def read_circular_table(path: str | Path) -> pd.DataFrame:
    """
    Read known ecDNA/circular structure tables.

    Supports both:
      1. Headered tables with AmpliconID/Chrom/Start/End/Strand columns.
      2. ASK2/ASK no-header circular TSV:
         Chrom, Start, End, Strand, Score1, Score2, AmpliconID,
         Gene, CancerGene, SegmentEdge.
    """
    df = read_table(path)
    required = {"AmpliconID", "Chrom", "Start", "End", "Strand"}
    if required.issubset(df.columns):
        return df

    raw = pd.read_csv(path, sep="\t", header=None)
    if raw.shape[1] < 7:
        raise ValueError(
            f"Cannot parse circular ecDNA structure table: {path}. "
            "Expected either a headered table with AmpliconID/Chrom/Start/End/Strand "
            "or an ASK no-header table with at least 7 columns."
        )
    names = ASK_CIRCULAR_NOHEADER_COLUMNS[: raw.shape[1]]
    if raw.shape[1] > len(names):
        names = names + [f"Extra{i}" for i in range(raw.shape[1] - len(names))]
    raw.columns = names
    return raw


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def write_pickle(obj, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(obj, handle)


def ask_style_output_paths(outdir: str | Path, bam_path: str | Path, outprefix: str | Path | None = None) -> dict[str, Path]:
    """
    Return ASK-compatible output paths for targeted ASK-search.

    The reference ASK runner writes files as <outprefix>_ask_*. This helper
    keeps the same structure while allowing --outdir to remain the simple CLI
    entry point for this module.
    """
    outdir = Path(outdir)
    prefix = Path(outprefix) if outprefix else outdir / Path(bam_path).stem
    return prepare_ask_output_paths(prefix)


def cleanup_legacy_bam_search_outputs(outdir: str | Path) -> None:
    """
    Remove old debug/auxiliary files from earlier ASK-search runs.
    """
    outdir = Path(outdir)
    patterns = [
        "target_*.tsv",
        "target_*.bedgraph",
        "target_*.txt",
        "segment_copy_number.tsv",
        "circle_copy_number_summary.tsv",
        "targeted_bppair_matches.tsv",
        "bam_junction_evidence.tsv",
        "breakpoint_support_summary.tsv",
        "*_ask_jcs_junctions.tsv",
        "cell_by_breakpoint_matrix.tsv",
        "cell_circle_summary.tsv",
        "cell_ecDNA_annotation.tsv",
    ]
    for pattern in patterns:
        for path in outdir.glob(pattern):
            if path.is_file():
                path.unlink()


def first_existing_path(paths: Iterable[str | Path | None]) -> Path | None:
    for path in paths:
        if path and Path(path).exists():
            return Path(path)
    return None


def prioritize_module_paths(paths: Iterable[str | Path | None]) -> None:
    """
    Put module search paths in the given priority order.
    """
    ordered: list[str] = []
    seen: set[str] = set()
    for path in paths:
        if not path:
            continue
        text = str(path)
        if Path(text).exists() and text not in seen:
            ordered.append(text)
            seen.add(text)
    for path in reversed(ordered):
        if path in sys.path:
            sys.path.remove(path)
        sys.path.insert(0, path)


def resolve_annotation_files(
    genome: str = "hg38",
    genefile: str | Path | None = None,
    sefile: str | Path | None = None,
    cgfile: str | Path | None = None,
) -> tuple[Path | None, Path | None, Path | None]:
    """
    Resolve ASK-style annotation files for Gene/CancerGene/SE annotation.
    """
    files = resolve_ask_annotation_files(
        genome=genome,
        genefile=genefile,
        sefile=sefile,
        cgfile=cgfile,
    )
    return files["genefile"], files["sefile"], files["cgfile"]


def annotate_amplicon_table(
    df: pd.DataFrame,
    genefile: str | Path | None = None,
    sefile: str | Path | None = None,
    cgfile: str | Path | None = None,
    ask_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    Add Gene, CancerGene, and SE columns with ASK grange/misc helpers.
    """
    if df.empty:
        out = df.copy()
        for col in ["Gene", "CancerGene", "SE"]:
            if col not in out.columns:
                out[col] = []
        return out
    modules = load_ask_modules(ask_path)
    import grange
    import misc

    out = df.copy()
    if genefile is not None:
        out = grange.GRange.map_gene(out, str(genefile))
    elif "Gene" not in out.columns:
        out["Gene"] = ""
    if cgfile is not None and "Gene" in out.columns:
        out = misc.map_cancergene(out, str(cgfile))
    elif "CancerGene" not in out.columns:
        out["CancerGene"] = ""
    if sefile is not None:
        out = grange.GRange.map_senchancer(out, str(sefile))
    elif "SE" not in out.columns:
        out["SE"] = ""
    return out


def write_junction_sequence_outputs(
    bp_pair: pd.DataFrame,
    bam_path: str | Path,
    output_align: str | Path,
    output_align_dir: str | Path,
    ask_path: str | Path | None = None,
    circ_anno: pd.DataFrame | None = None,
) -> None:
    """
    Write ASK-style breakpoint junction alignment/sequence outputs.
    """
    output_align = Path(output_align)
    output_align_dir = Path(output_align_dir)
    output_align.parent.mkdir(parents=True, exist_ok=True)
    output_align_dir.mkdir(parents=True, exist_ok=True)

    if bp_pair.empty:
        output_align.write_text("")
        return

    modules = load_ask_modules(ask_path)
    bpjoint = modules["bpjoint"]

    try:
        bpjoint.ouput_alignment(bp_pair, str(bam_path), str(output_align))
    except Exception:
        with output_align.open("w") as handle:
            for idx, row in bp_pair.reset_index(drop=True).iterrows():
                handle.write(
                    ">"
                    + "\t".join(
                        map(
                            str,
                            [
                                row["Chrom1"],
                                row["Coord1"],
                                row["Clip1"],
                                row["Chrom2"],
                                row["Coord2"],
                                row["Clip2"],
                                row.get("Count", ""),
                                row.get("offset", ""),
                            ],
                        )
                    )
                    + "\n"
                )
                handle.write(str(row.get("Seq", "")) + "\n\n")

    for old in output_align_dir.glob("junction_*.txt"):
        old.unlink()
    for old in output_align_dir.glob("circ_*.tsv"):
        old.unlink()

    if circ_anno is not None and not circ_anno.empty:
        try:
            bpjoint.ouput_alignment(bp_pair, str(bam_path), str(output_align_dir), circ_anno)
            if any(output_align_dir.glob("circ_*.tsv")):
                return
        except TypeError:
            pass
        except Exception:
            pass
        _write_circ_junction_sequence_outputs(bp_pair, str(bam_path), output_align_dir, circ_anno, bpjoint)
        return

    _write_breakpoint_junction_sequence_outputs(bp_pair, str(bam_path), output_align_dir, bpjoint)


def _write_breakpoint_junction_sequence_outputs(
    bp_pair: pd.DataFrame,
    bam_path: str,
    output_align_dir: Path,
    bpjoint: object,
) -> None:
    for idx, row in bp_pair.reset_index(drop=True).iterrows():
        stem = (
            f"junction_{idx:04d}_"
            f"{row['Chrom1']}_{int(row['Coord1'])}_{row['Clip1']}_"
            f"{row['Chrom2']}_{int(row['Coord2'])}_{row['Clip2']}"
        )
        _write_bp_alignment_file(output_align_dir / f"{stem}.txt", row, bam_path, bpjoint)


def _write_circ_junction_sequence_outputs(
    bp_pair: pd.DataFrame,
    bam_path: str,
    output_align_dir: Path,
    circ_anno: pd.DataFrame,
    bpjoint: object,
) -> None:
    bp_norm = bp_pair.reset_index(drop=True).copy()
    for amplicon_id, sub in circ_anno.groupby("AmpliconID", sort=False):
        matched = _bp_pairs_for_amplicon(bp_norm, sub)
        out = output_align_dir / f"{amplicon_id}.tsv"
        with out.open("w") as handle:
            for _, row in matched.iterrows():
                _write_bp_alignment_handle(handle, row, bam_path, bpjoint)


def _bp_pairs_for_amplicon(bp_pair: pd.DataFrame, circ_sub: pd.DataFrame) -> pd.DataFrame:
    circ = circ_sub.copy().reset_index(drop=True)
    if circ.empty:
        return bp_pair.iloc[0:0].copy()
    circ.loc[circ["Strand"].astype(str) == "-", ["Start", "End"]] = circ.loc[
        circ["Strand"].astype(str) == "-", ["End", "Start"]
    ].values
    keep: list[int] = []
    seen: set[int] = set()
    for idx in range(len(circ)):
        current = circ.iloc[idx]
        nxt = circ.iloc[(idx + 1) % len(circ)]
        expected = [
            (current["Chrom"], int(current["End"]), "R", nxt["Chrom"], int(nxt["Start"]), "L"),
        ]
        if len(circ) == 1:
            expected.append((current["Chrom"], int(current["Start"]), "L", current["Chrom"], int(current["End"]), "R"))
        for row_idx, row in bp_pair.iterrows():
            if row_idx in seen:
                continue
            if any(_same_bp_pair(row, pair) for pair in expected):
                keep.append(row_idx)
                seen.add(row_idx)
    return bp_pair.loc[keep].copy()


def _same_bp_pair(row: pd.Series, pair: tuple[object, int, str, object, int, str], window: int = 2) -> bool:
    a = (str(pair[0]), int(pair[1]), str(pair[2]))
    b = (str(pair[3]), int(pair[4]), str(pair[5]))
    row_a = (str(row["Chrom1"]), int(row["Coord1"]), str(row["Clip1"]))
    row_b = (str(row["Chrom2"]), int(row["Coord2"]), str(row["Clip2"]))
    return (_same_breakend(row_a, a, window) and _same_breakend(row_b, b, window)) or (
        _same_breakend(row_a, b, window) and _same_breakend(row_b, a, window)
    )


def _same_breakend(a: tuple[str, int, str], b: tuple[str, int, str], window: int) -> bool:
    return a[0] == b[0] and a[2] == b[2] and abs(a[1] - b[1]) <= window


def _write_bp_alignment_file(out: Path, row: pd.Series, bam_path: str, bpjoint: object) -> None:
    with out.open("w") as handle:
        _write_bp_alignment_handle(handle, row, bam_path, bpjoint)


def _write_bp_alignment_handle(handle: object, row: pd.Series, bam_path: str, bpjoint: object) -> None:
    handle.write(
        ">"
        + "\t".join(
            map(
                str,
                [
                    row["Chrom1"],
                    row["Coord1"],
                    row["Clip1"],
                    row["Chrom2"],
                    row["Coord2"],
                    row["Clip2"],
                    row.get("Count", ""),
                    row.get("offset", ""),
                ],
            )
        )
        + "\n"
    )
    if str(row.get("Seq", "")) == "PE_Support":
        handle.write("PE_Support\n\n")
        return
    try:
        alignment = bpjoint.get_alignment(
            bam_path,
            (row["Chrom1"], int(row["Coord1"]), row["Clip1"]),
            (row["Chrom2"], int(row["Coord2"]), row["Clip2"]),
            offset=int(row.get("offset", 0)),
        )
        for line in alignment:
            handle.write(f"{line}\n")
        handle.write("\n\n")
    except Exception:
        handle.write(str(row.get("Seq", "")) + "\n\n")


def simple_amplicon_plot(
    amp_df: pd.DataFrame,
    bin_norm: pd.DataFrame,
    fig_dir: str | Path,
    binsize: int = 10000,
    max_plots: int = 20,
) -> None:
    """
    Lightweight fallback plot when ASK trackplot cannot render.
    """
    if amp_df.empty:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    for i, (amplicon_id, sub) in enumerate(amp_df.groupby("AmpliconID", sort=False)):
        if i >= max_plots:
            break
        chroms = sorted(set(sub["Chrom"].astype(str)))
        start = int(sub["Start"].min())
        end = int(sub["End"].max())
        fig, ax = plt.subplots(figsize=(12, 3.5))
        for _, seg in sub.iterrows():
            ax.plot([seg["Start"], seg["End"]], [1, 1], linewidth=8, solid_capstyle="butt")
            ax.text((seg["Start"] + seg["End"]) / 2, 1.08, str(seg.get("Gene", ""))[:30], ha="center", fontsize=7)
        for j in range(len(sub)):
            a = sub.iloc[j]
            b = sub.iloc[(j + 1) % len(sub)]
            x1 = int(a["End"])
            x2 = int(b["Start"])
            center = (x1 + x2) / 2
            width = max(abs(x2 - x1), binsize)
            ax.add_patch(Arc((center, 1.02), width=width, height=0.55, theta1=0, theta2=180, linewidth=1.2))
        if not bin_norm.empty:
            cn = bin_norm[
                (bin_norm["Chrom"].astype(str).isin(chroms))
                & (bin_norm["Coord"].astype(int) >= start - 5 * binsize)
                & (bin_norm["Coord"].astype(int) <= end + 5 * binsize)
            ]
            if not cn.empty:
                ax2 = ax.twinx()
                ax2.plot(cn["Coord"], cn["CN"], color="0.35", linewidth=0.8, alpha=0.8)
                ax2.set_ylabel("CN")
        ax.set_title(str(amplicon_id))
        ax.set_xlabel(";".join(chroms))
        ax.set_yticks([])
        ax.set_xlim(start - max(binsize, int((end - start) * 0.05)), end + max(binsize, int((end - start) * 0.05)))
        fig.tight_layout()
        fig.savefig(fig_dir / f"circular_{amplicon_id}.pdf")
        plt.close(fig)


def plot_target_amplicons(
    circ_anno: pd.DataFrame,
    line_anno: pd.DataFrame,
    cn_amp: pd.DataFrame,
    genefile: str | Path | None,
    sefile: str | Path | None,
    bin_norm: pd.DataFrame,
    binsize: int,
    fig_dir: str | Path,
    ask_path: str | Path | None = None,
) -> None:
    """
    Plot amplicon figures with the original ASK plotting function.
    """
    modules = load_ask_modules(ask_path)
    ask_core = modules["ask"]
    import trackplot

    print(f"plot amplicons - ask.py: {ask_core.__file__}")
    print(f"plot amplicons - trackplot.py: {trackplot.__file__}")
    print(f"plot amplicons - genefile: {genefile}")
    plot_ask_amplicons(
        ask_core,
        circ_anno,
        line_anno,
        cn_amp,
        genefile,
        sefile,
        bin_norm,
        fig_dir=fig_dir,
        binsize=binsize,
        plot_n=5,
        fig_width=15,
        fontsize=12,
        ext=0.3,
        clean_old_pdf=True,
    )


def parse_sample_arg(value: str) -> SampleInput:
    """
    Parse one sample argument.

    Accepted forms:
      sample_id=/path/to/circular.tsv
      sample_id,filename=/path/to/circular.tsv

    The second form is useful when the ASK circular file should keep the
    original EGA filename used by downstream Seurat annotation code.
    """
    if "=" not in value:
        raise ValueError(
            "--circular must be sample_id=path or sample_id,filename=path"
        )
    left, right = value.split("=", 1)
    if "," in left:
        sample_id, filename = left.split(",", 1)
    else:
        sample_id, filename = left, left
    return SampleInput(sample_id=sample_id, filename=filename, circular_path=Path(right))


def load_circular_tables(samples: Iterable[SampleInput]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for sample in samples:
        if isinstance(sample.circular_path, pd.DataFrame):
            df = sample.circular_path.copy()
        else:
            df = read_circular_table(sample.circular_path)
        df["SampleID"] = sample.sample_id
        df["Filename"] = sample.filename or sample.sample_id
        frames.append(df)
    if not frames:
        raise ValueError("No circular amplicon files were provided.")
    out = pd.concat(frames, ignore_index=True)
    require_columns(out, ["AmpliconID", "Chrom", "Start", "End", "Strand", "Filename"])
    out["CircleID"] = out["AmpliconID"].astype(str) + "_" + out["Filename"].astype(str)
    return out


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(";") if item.strip()]


def summarize_cancer_genes(values: pd.Series, target_genes: Iterable[str] | None = None) -> str:
    genes: list[str] = []
    for value in values.dropna():
        genes.extend(split_genes(value))
    if target_genes:
        targets = set(target_genes)
        genes = [gene for gene in genes if gene in targets]
    return "; ".join(sorted(set(genes)))


def extract_breakpoint_pairs(circular_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert ASK circular amplicon segments into circular breakpoint pairs.

    For segments on the negative strand, Start/End are swapped before pairing.
    Within each CircleID, each segment end is connected to the next segment start;
    the last segment closes the circle by connecting back to the first segment.
    """
    require_columns(circular_df, ["CircleID", "AmpliconID", "Filename", "Chrom", "Start", "End", "Strand"])
    df = circular_df.copy()
    minus = df["Strand"].astype(str) == "-"
    df.loc[minus, ["Start", "End"]] = df.loc[minus, ["End", "Start"]].to_numpy()

    rows: list[dict[str, object]] = []
    for circle_id, circle_df in df.groupby("CircleID", sort=False):
        circle_df = circle_df.reset_index(drop=True)
        n_segments = len(circle_df)
        cancer_gene = summarize_cancer_genes(circle_df.get("CancerGene", pd.Series(dtype=object)))
        all_gene = summarize_cancer_genes(circle_df.get("Gene", pd.Series(dtype=object)))

        for i, row in circle_df.iterrows():
            next_row = circle_df.iloc[(i + 1) % n_segments]
            bp_id = (
                f"{circle_id}|{row['Chrom']}:{int(row['End'])}"
                f"->{next_row['Chrom']}:{int(next_row['Start'])}"
            )
            rows.append(
                {
                    "CircleID": circle_id,
                    "AmpliconID": row["AmpliconID"],
                    "Filename": row["Filename"],
                    "SampleID": row.get("SampleID", row["Filename"]),
                    "BreakpointID": bp_id,
                    "Chrom1": row["Chrom"],
                    "Coord1": int(row["End"]),
                    "Chrom2": next_row["Chrom"],
                    "Coord2": int(next_row["Start"]),
                    "Strand1": row.get("Strand", ""),
                    "Strand2": next_row.get("Strand", ""),
                    "Gene": all_gene,
                    "CancerGene": cancer_gene,
                    "SegmentCount": n_segments,
                }
            )

    return pd.DataFrame(rows).drop_duplicates(
        ["CircleID", "Chrom1", "Coord1", "Chrom2", "Coord2"]
    )


def infer_clip_from_known_breakpoints(breakpoints: pd.DataFrame) -> pd.DataFrame:
    """
    Add ASK-style Clip1/Clip2 labels for known ecDNA junctions.

    A circular segment contributes its right boundary (R) as Coord1 and the next
    segment contributes its left boundary (L) as Coord2. This matches the ASK
    convention used by getbppair/bppair.
    """
    out = breakpoints.copy()
    if "Clip1" not in out.columns:
        out["Clip1"] = "R"
    if "Clip2" not in out.columns:
        out["Clip2"] = "L"
    return out


def canonical_breakpoint_key(chrom1: object, coord1: object, chrom2: object, coord2: object) -> tuple:
    a = (str(chrom1), int(coord1))
    b = (str(chrom2), int(coord2))
    return tuple(sorted((a, b)))


def add_breakpoint_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_bp_key"] = [
        canonical_breakpoint_key(row.Chrom1, row.Coord1, row.Chrom2, row.Coord2)
        for row in out.itertuples()
    ]
    return out


def parse_list_cell(value: object) -> list:
    if isinstance(value, (list, tuple, set)):
        return list(value)
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return value.tolist()
    if pd.isna(value):
        return []
    text = str(value)
    for parser in (json.loads, ast.literal_eval):
        try:
            parsed = parser(text)
            if isinstance(parsed, list):
                return parsed
        except Exception:
            pass
    return [item.strip() for item in text.strip("[]").split(",") if item.strip()]


def safe_get_tag(read, tag: str, default=None):
    try:
        return read.get_tag(tag)
    except Exception:
        return default


def ask_check_read(read, mapq: int = 20, nmmax: int = 1) -> bool:
    """
    ASK-compatible read filter.

    ASK uses mapped, MAPQ-filtered, low-NM, non-duplicate reads for breakpoint
    evidence. Reads without NM are allowed if they otherwise pass filters, which
    makes this wrapper usable across aligners that omit NM.
    """
    if read.is_unmapped or read.is_duplicate or read.mapping_quality < mapq:
        return False
    nm = safe_get_tag(read, "NM", 0)
    return nm is None or nm <= nmmax


def read_clip_side(read) -> tuple[str, int] | tuple[None, None]:
    """
    Return ASK-style clip side and coordinate for one read.

    L: left soft/hard clipped read at reference_start.
    R: right soft/hard clipped read at reference_end - 1.
    """
    if read.cigartuples is None:
        return None, None
    if read.cigartuples[0][0] in {4, 5}:
        return "L", read.reference_start
    if read.cigartuples[-1][0] in {4, 5}:
        return "R", read.reference_end - 1
    return None, None


def read_barcode(read, tags: Iterable[str] = ("CB", "CR", "BX", "BC")):
    for tag in tags:
        value = safe_get_tag(read, tag)
        if value is not None:
            return value
    return None


def parse_sa_tag(read) -> list[dict[str, object]]:
    """
    Parse SA tag entries into dictionaries.

    SA format: rname,pos,strand,CIGAR,mapQ,NM; pos is 1-based.
    """
    sa = safe_get_tag(read, "SA")
    if not sa:
        return []
    rows = []
    for item in sa.split(";"):
        if not item:
            continue
        fields = item.split(",")
        if len(fields) < 6:
            continue
        rows.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]) - 1,
                "strand": fields[2],
                "cigar": fields[3],
                "mapq": int(fields[4]),
                "nm": int(fields[5]),
            }
        )
    return rows


def load_ask_getbppair(ask_path: str | Path | None = None):
    """
    Load ASK getbppair module without requiring it to be installed.

    getbppair.get_bp_pair is the ASK step that converts SA/supplementary
    alignments into breakpoint-pair candidates and returns the junction
    sequence in the Seq column.
    """
    candidate_paths = []
    if ask_path:
        candidate_paths.append(str(ask_path))
    candidate_paths.extend(DEFAULT_ASK_PATHS)
    prioritize_module_paths(candidate_paths)
    import getbppair

    return getbppair


def load_ask_modules(ask_path: str | Path | None = None):
    """
    Load ASK modules used by the targeted pipeline.
    """
    candidate_paths = []
    if ask_path:
        candidate_paths.append(str(ask_path))
    candidate_paths.extend(DEFAULT_ASK_PATHS)
    prioritize_module_paths(candidate_paths)
    import bpdetect
    import bpjoint
    import bppair
    import cbs
    import cndetect
    import ggraph
    import getbppair

    ask_file = first_existing_path([Path(path) / "ask.py" for path in candidate_paths])
    if ask_file is None:
        raise FileNotFoundError("Could not find ASK core module ask.py")
    spec = importlib.util.spec_from_file_location("ask_core_module", ask_file)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load ASK core module from {ask_file}")
    ask_core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ask_core)

    return {
        "ask": ask_core,
        "bpdetect": bpdetect,
        "bpjoint": bpjoint,
        "bppair": bppair,
        "cbs": cbs,
        "cndetect": cndetect,
        "ggraph": ggraph,
        "getbppair": getbppair,
    }


def list2df_count_simple(items: list[tuple], columns: list[str]) -> pd.DataFrame:
    if not items:
        return pd.DataFrame(columns=columns + ["Count"])
    return (
        pd.DataFrame(items, columns=columns)
        .groupby(columns)
        .size()
        .reset_index(name="Count")
        .sort_values("Count", ascending=False)
        .reset_index(drop=True)
    )


def target_chromosomes_from_circular(circular_df: pd.DataFrame) -> list[str]:
    return sorted(circular_df["Chrom"].astype(str).unique())


def known_single_segment_circles(circular_df: pd.DataFrame) -> pd.DataFrame:
    if circular_df.empty or "CircleID" not in circular_df.columns:
        return pd.DataFrame(columns=["Chrom", "Start", "End"])
    counts = circular_df.groupby("CircleID")["CircleID"].transform("size")
    return circular_df.loc[counts == 1, ["Chrom", "Start", "End"]].drop_duplicates().copy()


def target_process_alignment(
    bam_path: str | Path,
    breakpoints: pd.DataFrame,
    target_chroms: Iterable[str],
    circular_df: pd.DataFrame | None = None,
    cn_background: str = "genome",
    window: int = 200,
    binsize: int = 10000,
    mapq: int = 20,
    nmmax: int = 3,
    ask_path: str | Path | None = None,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """
    Targeted replacement for ASK process_alignment.

    Only breakpoint windows are scanned for clip/SA evidence; only chromosomes
    containing the known ecDNA structure are scanned for bin counts.
    """
    modules = load_ask_modules(ask_path)
    getbppair = modules["getbppair"]

    left_clip: list[tuple] = []
    right_clip: list[tuple] = []
    reads_with_sa: dict = {}
    seen = set()
    bps = infer_clip_from_known_breakpoints(breakpoints)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for bp in bps.itertuples():
            for read in fetch_reads_around_breakpoint(
                bam, bp.Chrom1, bp.Coord1, bp.Chrom2, bp.Coord2, window
            ):
                key = (read.query_name, read.flag, read.reference_id, read.reference_start)
                if key in seen:
                    continue
                seen.add(key)
                add_read_to_ask_sa_group(read, reads_with_sa, mapq=0, nmmax=nmmax)
                if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
                    continue
                clip, coord = read_clip_side(read)
                if clip == "L":
                    left_clip.append((read.reference_name, int(coord), "L"))
                elif clip == "R":
                    right_clip.append((read.reference_name, int(coord), "R"))

    bp_all = {
        "left": list2df_count_simple(left_clip, ["Chrom", "Coord", "Clip"]),
        "right": list2df_count_simple(right_clip, ["Chrom", "Coord", "Clip"]),
    }
    clip_bg = clip_to_bedgraph(bp_all)
    bin_count = target_region_count(
        bam_path,
        target_chroms=target_chroms,
        exclude_intervals=circular_df,
        cn_background=cn_background,
        binsize=binsize,
        mapq=mapq,
        nmmax=nmmax,
    )
    bp_duo_bam = getbppair.get_bp_pair(reads_with_sa, mapq=mapq) if reads_with_sa else pd.DataFrame()
    return bp_all, clip_bg, bin_count, bp_duo_bam, reads_with_sa


def clip_to_bedgraph(bp_all: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for side, sign in [("left", -1), ("right", 1)]:
        df = bp_all.get(side, pd.DataFrame()).copy()
        if df.empty:
            continue
        for row in df.itertuples():
            rows.append([row.Chrom, int(row.Coord), int(row.Coord) + 1, sign * int(row.Count)])
    return pd.DataFrame(rows, columns=["Chrom", "Start", "End", "Count"]).sort_values(["Chrom", "Start", "End"])


def target_region_count(
    bam_path: str | Path,
    target_chroms: Iterable[str],
    exclude_intervals: pd.DataFrame | None = None,
    cn_background: str = "genome",
    binsize: int = 10000,
    mapq: int = 20,
    nmmax: int = 3,
) -> pd.DataFrame:
    counts: dict[tuple[str, int], int] = {}
    norm_counts: dict[tuple[str, int], int] = {}
    target = set(map(str, target_chroms))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if cn_background == "genome":
            iterator = bam.fetch()
        else:
            iterators = []
            for chrom in target:
                try:
                    iterators.append(bam.fetch(chrom))
                except ValueError:
                    continue
            iterator = itertools.chain.from_iterable(iterators)
        for read in iterator:
            if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
                continue
            coord = int(read.reference_start // binsize * binsize)
            key = (read.reference_name, coord)
            norm_counts[key] = norm_counts.get(key, 0) + 1
            if str(read.reference_name) in target:
                counts[key] = counts.get(key, 0) + 1
    if not counts:
        return pd.DataFrame(columns=["Chrom", "Coord", "Count", "CN"])
    df = pd.DataFrame(
        [(chrom, coord, count) for (chrom, coord), count in counts.items()],
        columns=["Chrom", "Coord", "Count"],
    ).sort_values(["Chrom", "Coord"])
    norm_source = (
        pd.DataFrame(
            [(chrom, coord, count) for (chrom, coord), count in norm_counts.items()],
            columns=["Chrom", "Coord", "Count"],
        )
        if norm_counts
        else df
    )
    if cn_background != "genome" and exclude_intervals is not None and not exclude_intervals.empty:
        intervals = exclude_intervals[["Chrom", "Start", "End"]].drop_duplicates().copy()
        intervals["Chrom"] = intervals["Chrom"].astype(str)
        intervals["Start"] = intervals["Start"].astype(int)
        intervals["End"] = intervals["End"].astype(int)
        in_excluded = pd.Series(False, index=df.index)
        for chrom, sub in intervals.groupby("Chrom", sort=False):
            chrom_mask = df["Chrom"].astype(str) == chrom
            if not chrom_mask.any():
                continue
            for row in sub.itertuples(index=False):
                in_excluded |= chrom_mask & (df["Coord"].astype(int) < int(row.End)) & (
                    df["Coord"].astype(int) + binsize > int(row.Start)
                )
        if (~in_excluded).sum() >= 10:
            norm_source = df.loc[~in_excluded]
    mean_count = norm_source["Count"].mean()
    df["CN"] = 2 * df["Count"] / mean_count if mean_count else 0
    return df.reset_index(drop=True)


def target_detect_amplified_segment(
    bin_count: pd.DataFrame,
    bp_all: dict[str, pd.DataFrame],
    circular_df: pd.DataFrame,
    binsize: int = 10000,
    ask_path: str | Path | None = None,
    min_cn: float | None = None,
    min_segsize_scale: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Targeted replacement for ASK detect_amplified_segment on ecDNA chromosomes.

    If CBS/cndetect does not return an amplified segment, the known ecDNA
    segment intervals are used as a fallback cn_amp so downstream DFS can still
    evaluate junction support.
    """
    modules = load_ask_modules(ask_path)
    cndetect = modules["cndetect"]
    cbs = modules["cbs"]

    if bin_count.empty:
        cn_amp = known_segments_to_cn_amp(circular_df, bin_count, bp_all, binsize)
        return cn_amp, pd.DataFrame(), bin_count

    bin_norm = cndetect.smooth_cn(bin_count.copy())
    try:
        cn_seg = cbs.cbs(bin_norm, binsize=binsize)
        cn_amp_raw = cbs.cbs_amplicon(cn_seg, min_cn=min_cn)
        cn_amp_raw = cndetect.segment_filter_by_size(
            cn_amp_raw, binsize=binsize, fold=min_segsize_scale
        )
        cn_amp = cndetect.finemap_cn_segment_boundary(cn_amp_raw, bp_all, binsize=binsize)
    except Exception:
        cn_seg = pd.DataFrame()
        cn_amp = pd.DataFrame()

    if cn_amp.empty:
        cn_amp = known_segments_to_cn_amp(circular_df, bin_norm, bp_all, binsize)
    else:
        cn_amp = ensure_cn_amp_annotation(cn_amp, circular_df)
    return cn_amp, cn_seg, bin_norm


def known_segments_to_cn_amp(
    circular_df: pd.DataFrame,
    bin_norm: pd.DataFrame,
    bp_all: dict[str, pd.DataFrame],
    binsize: int = 10000,
) -> pd.DataFrame:
    rows = []
    grouped = circular_df.groupby(["Chrom", "Start", "End"], sort=False)
    for (chrom, start, end), sub in grouped:
        cn, _ = weighted_cn_for_interval(bin_norm, chrom, start, end, binsize)
        clip_left = clip_count_at(bp_all, chrom, start, "L")
        clip_right = clip_count_at(bp_all, chrom, end, "R")
        rows.append(
            {
                "Chrom": chrom,
                "Start": int(start),
                "End": int(end),
                "Count": 0,
                "CN": 2 if cn is None else cn,
                "ClipLeft": clip_left,
                "ClipRight": clip_right,
                "Gene": summarize_cancer_genes(sub.get("Gene", pd.Series(dtype=object))),
                "CancerGene": summarize_cancer_genes(sub.get("CancerGene", pd.Series(dtype=object))),
            }
        )
    return pd.DataFrame(rows)


def clip_count_at(bp_all: dict[str, pd.DataFrame], chrom: object, coord: object, clip: str, window: int = 2) -> int:
    side = "left" if clip == "L" else "right"
    df = bp_all.get(side, pd.DataFrame())
    if df.empty:
        return 0
    sub = df[
        (df["Chrom"].astype(str) == str(chrom))
        & (df["Clip"].astype(str) == clip)
        & ((df["Coord"].astype(int) - int(coord)).abs() <= window)
    ]
    return int(sub["Count"].sum()) if not sub.empty else 0


def ensure_cn_amp_annotation(cn_amp: pd.DataFrame, circular_df: pd.DataFrame) -> pd.DataFrame:
    out = cn_amp.copy()
    for col in ["Gene", "CancerGene"]:
        if col not in out.columns:
            out[col] = ""
    return out


def known_breakpoint_seed_table(breakpoints: pd.DataFrame) -> pd.DataFrame:
    """
    Build ASK bp_duo_bam-like seeds from known ecDNA breakpoint pairs.

    bpjoint.get_bppair reads this table by column position, so the order is
    intentionally fixed to match ASK's bp_duo_bam layout.
    """
    bps = infer_clip_from_known_breakpoints(breakpoints)
    require_columns(bps, ["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2"])
    rows = []
    for row in bps.itertuples():
        rows.append(
            {
                "Chrom1": row.Chrom1,
                "Coord1": int(row.Coord1),
                "Clip1": row.Clip1,
                "Chrom2": row.Chrom2,
                "Coord2": int(row.Coord2),
                "Clip2": row.Clip2,
                "Count": 1,
                "offset": 0,
                "Seq": "Known_Input",
                "Readsname": [],
                "CircleID": getattr(row, "CircleID", ""),
                "BreakpointID": getattr(row, "BreakpointID", ""),
            }
        )
    columns = [
        "Chrom1",
        "Coord1",
        "Clip1",
        "Chrom2",
        "Coord2",
        "Clip2",
        "Count",
        "offset",
        "Seq",
        "Readsname",
        "CircleID",
        "BreakpointID",
    ]
    return pd.DataFrame(rows, columns=columns)


def bpduo_core_columns(df: pd.DataFrame) -> pd.DataFrame:
    columns = ["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", "Count", "offset", "Seq"]
    if df.empty:
        return pd.DataFrame(columns=columns)
    out = df.copy()
    for col, default in [("Count", 1), ("offset", 0), ("Seq", "")]:
        if col not in out.columns:
            out[col] = default
    return out[columns]


def bpduo_bam_seed_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep columns and order expected by bpjoint.get_bppair(type=2).
    """
    columns = [
        "Chrom1",
        "Coord1",
        "Clip1",
        "Chrom2",
        "Coord2",
        "Clip2",
        "Count",
        "offset",
        "Seq",
        "Readsname",
    ]
    if df.empty:
        return pd.DataFrame(columns=columns)
    out = df.copy()
    for col, default in [("Count", 1), ("offset", 0), ("Seq", ""), ("Readsname", [])]:
        if col not in out.columns:
            out[col] = default
    return out[columns]


def breakpoint_pair_close(row: pd.Series, seed: pd.Series, window: int) -> bool:
    direct = (
        breakpoint_end_matches(row["Chrom1"], row["Coord1"], row["Clip1"], seed["Chrom1"], seed["Coord1"], seed["Clip1"], window)
        and breakpoint_end_matches(row["Chrom2"], row["Coord2"], row["Clip2"], seed["Chrom2"], seed["Coord2"], seed["Clip2"], window)
    )
    reverse = (
        breakpoint_end_matches(row["Chrom1"], row["Coord1"], row["Clip1"], seed["Chrom2"], seed["Coord2"], seed["Clip2"], window)
        and breakpoint_end_matches(row["Chrom2"], row["Coord2"], row["Clip2"], seed["Chrom1"], seed["Coord1"], seed["Clip1"], window)
    )
    return bool(direct or reverse)


def filter_bp_duo_to_known_windows(bp_duo: pd.DataFrame, known_seed: pd.DataFrame, window: int = 200) -> pd.DataFrame:
    if bp_duo.empty or known_seed.empty:
        return pd.DataFrame(columns=bp_duo.columns)
    rows = []
    for _, cand in bp_duo.iterrows():
        if any(breakpoint_pair_close(cand, seed, window) for _, seed in known_seed.iterrows()):
            rows.append(cand)
    subset = ["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2"]
    return pd.DataFrame(rows).drop_duplicates(subset=subset).reset_index(drop=True) if rows else pd.DataFrame(columns=bp_duo.columns)


def target_detect_bp_pair(
    bam_path: str | Path,
    known_seed: pd.DataFrame,
    bp_duo_bam: pd.DataFrame,
    bp_all: dict[str, pd.DataFrame],
    cn_amp: pd.DataFrame,
    binsize: int = 10000,
    window: int = 200,
    min_nt: int = 5,
    mapq: int = 20,
    nmmax: int = 3,
    ask_path: str | Path | None = None,
    include_sa_candidates: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Targeted ASK breakpoint-pair step.

    Known breakpoint pairs are supplied as bp_duo_bam seeds, so bpjoint only
    performs local matching/junction sequence recovery for those pairs. SA
    candidates found in the breakpoint windows can optionally be retained if
    they fall near the known pairs.
    """
    modules = load_ask_modules(ask_path)
    bpjoint = modules["bpjoint"]
    bpdetect = modules["bpdetect"]

    seed = bpduo_bam_seed_columns(known_seed)
    assembled_pairs = []
    empty_seed = pd.DataFrame(columns=["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", "Count", "offset", "Seq", "Readsname"])
    for _, seed_row in seed.iterrows():
        bp_cand_df = known_pair_to_bp_cand_df(seed_row)
        try:
            pair_df = bpjoint.get_bppair(
                str(bam_path),
                bp_cand_df,
                empty_seed,
                min_nt=min_nt,
                multi_process=False,
                num_cores=1,
            )
        except Exception:
            pair_df = pd.DataFrame()
        if not pair_df.empty:
            assembled_pairs.append(pair_df)

    if include_sa_candidates and not bp_duo_bam.empty:
        sa_near_known = filter_bp_duo_to_known_windows(bpduo_bam_seed_columns(bp_duo_bam), seed, window=window)
        if not sa_near_known.empty:
            sa_pairs = bpjoint.get_bppair(
                str(bam_path),
                pd.DataFrame(columns=["Chrom", "Coord", "Clip", "Group"]),
                sa_near_known,
                min_nt=min_nt,
                multi_process=False,
                num_cores=1,
            )
            if not sa_pairs.empty:
                assembled_pairs.append(sa_pairs)
    else:
        sa_near_known = pd.DataFrame()

    if assembled_pairs:
        bp_duo = pd.concat([bpduo_core_columns(df) for df in assembled_pairs], ignore_index=True)
        bp_duo = (
            bp_duo.sort_values("Count", ascending=False)
            .drop_duplicates(["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2"])
            .reset_index(drop=True)
        )
    else:
        bp_duo = pd.DataFrame(columns=["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", "Count", "offset", "Seq"])

    if bp_duo.empty:
        bp_cand_stats = seed_to_bp_cand_stats(seed, bp_all)
    else:
        bp_cand_stats = bpdetect.refine_bp_cand_bam(bp_duo, bp_all)
    return bp_duo, bp_cand_stats, sa_near_known


def seed_to_bp_cand_stats(seed: pd.DataFrame, bp_all: dict[str, pd.DataFrame]) -> pd.DataFrame:
    cols = ["Chrom", "Coord", "Clip"]
    if seed.empty:
        return pd.DataFrame(columns=cols + ["CleanBP", "ClipDepth", "InDepth", "OutDepth"])
    bp_df = pd.concat(
        [
            seed[["Chrom1", "Coord1", "Clip1"]].set_axis(cols, axis=1),
            seed[["Chrom2", "Coord2", "Clip2"]].set_axis(cols, axis=1),
        ],
        ignore_index=True,
    ).drop_duplicates()
    rows = []
    for row in bp_df.itertuples():
        rows.append(
            {
                "Chrom": row.Chrom,
                "Coord": int(row.Coord),
                "Clip": row.Clip,
                "CleanBP": True,
                "ClipDepth": clip_count_at(bp_all, row.Chrom, row.Coord, row.Clip, window=2),
                "InDepth": 0,
                "OutDepth": 0,
            }
        )
    return pd.DataFrame(rows)


def known_pair_to_bp_cand_df(seed_row: pd.Series) -> pd.DataFrame:
    """
    Convert one known breakpoint pair into the bp_cand_df shape ASK assembles.
    """
    return pd.DataFrame(
        [
            {
                "Chrom": seed_row["Chrom1"],
                "Coord": int(seed_row["Coord1"]),
                "Clip": seed_row["Clip1"],
                "Group": 1,
            },
            {
                "Chrom": seed_row["Chrom2"],
                "Coord": int(seed_row["Coord2"]),
                "Clip": seed_row["Clip2"],
                "Group": 1,
            },
        ]
    )


def empty_amplicon_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    circ_cols = ["Chrom", "Start", "End", "Strand", "SplitCount", "CN", "AmpliconID", "Gene"]
    score_cols = [
        "AmpliconID",
        "Chrom1",
        "Start",
        "End",
        "Seg_num",
        "Length",
        "SplitCount_sum",
        "SplitCount_mean",
        "SplitCount_std",
        "CN_sum",
        "CN_mean",
        "CN_std",
        "Left_CN",
        "Right_CN",
        "Gene_num",
        "Cancergene_num",
        "SE_num",
        "F1",
        "F2",
        "F3",
        "F4",
        "Score",
    ]
    return (
        pd.DataFrame(columns=circ_cols),
        pd.DataFrame(columns=circ_cols),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(columns=score_cols),
        pd.DataFrame(),
    )


def target_construct_amplicon(
    bp_duo: pd.DataFrame,
    bp_cand_stats: pd.DataFrame,
    cn_amp: pd.DataFrame,
    bin_norm: pd.DataFrame,
    known_single_segments: pd.DataFrame | None = None,
    min_junc_cnt: int = 1,
    bpp_min_dist: int = 50,
    segment_restrict: bool = False,
    subseg: bool = True,
    rm_amp_df: bool = False,
    ask_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    ASK-style DFS construction from targeted breakpoint pairs.
    """
    if bp_duo.empty:
        return empty_amplicon_outputs()

    modules = load_ask_modules(ask_path)
    bppair = modules["bppair"]
    ggraph = modules["ggraph"]

    cn_cols = ["Chrom", "Start", "End", "CN", "ClipLeft", "ClipRight", "Gene", "CancerGene"]
    cn_amp_use = cn_amp.copy()
    for col in cn_cols:
        if col not in cn_amp_use.columns:
            cn_amp_use[col] = "" if col in {"Gene", "CancerGene"} else 0
    cn_amp_use = cn_amp_use[cn_cols]

    bp_pair = bppair.bp_pair_filter(bpduo_core_columns(bp_duo), bpp_min_dist, min_junc_cnt)
    if bp_pair.empty:
        empty = empty_amplicon_outputs()
        return empty[0], empty[1], bp_pair, empty[3], empty[4], empty[5]

    bp_fine, cn_amp_merged = bppair.bp_refine(bp_pair, bp_cand_stats, cn_amp_use)
    seg = bppair.get_segment(
        bp_fine,
        bp_pair,
        cn_amp_merged,
        restrict=segment_restrict,
        subseg=subseg,
        rm_amp_df=rm_amp_df,
    )
    if not seg.empty:
        seg = bppair.add_cn(seg, cn_amp_use, bin_norm)

    gg = ggraph.GGraph()
    gg.build_ggraph_from_bp(bp_pair, bp_fine, seg)
    all_circ_path = gg.get_all_path_contain_circle()
    circ_uniq = gg.path_unique(all_circ_path)
    circ = gg.get_representitive_path(circ_uniq)
    circ_anno = gg.make_amplicon_df(circ)
    if getattr(gg, "single_segment_loop", None) and known_single_segments is not None and not known_single_segments.empty:
        col_sel = ["Chrom", "Start", "End", "Strand", "SplitCount", "CN"]
        single_circ = pd.DataFrame(gg.single_segment_loop, columns=col_sel)
        if not single_circ.empty:
            single_key = set(
                zip(
                    known_single_segments["Chrom"].astype(str),
                    known_single_segments["Start"].astype(int),
                    known_single_segments["End"].astype(int),
                )
            )
            single_circ = single_circ[
                [
                    (str(row.Chrom), int(row.Start), int(row.End)) in single_key
                    for row in single_circ.itertuples()
                ]
            ].copy()
        if not single_circ.empty:
            if not circ_anno.empty:
                existing = circ_anno[col_sel].drop_duplicates()
                single_circ = single_circ.merge(existing, on=col_sel, how="left", indicator=True)
                single_circ = single_circ[single_circ["_merge"] == "left_only"].drop(columns="_merge")
                next_index = circ_anno["AmpliconID"].nunique()
            else:
                next_index = 0
            if not single_circ.empty:
                single_circ = single_circ.reset_index(drop=True)
                single_circ["AmpliconID"] = [f"circ_{next_index + i}" for i in range(single_circ.shape[0])]
                circ_anno = pd.concat([circ_anno, single_circ], ignore_index=True)

    line_uniq = gg.path_unique(all_circ_path, type="linear")
    line = gg.get_representitive_path(line_uniq, type="linear")
    line_anno = gg.make_amplicon_df(line, type="linear")

    if not circ_anno.empty and not bin_norm.empty:
        try:
            circ_score = ggraph.add_stats_circ(circ_anno, bin_norm, binsize=10000)
        except Exception:
            circ_score = pd.DataFrame()
    else:
        circ_score = pd.DataFrame()

    return circ_anno, line_anno, bp_pair, seg, circ_score, bp_fine


def add_read_to_ask_sa_group(read, reads_with_sa: dict, mapq: int = 0, nmmax: int = 3) -> None:
    """
    Build the reads_with_SA structure expected by ASK getbppair.get_bp_pair.
    """
    if not read.has_tag("SA"):
        return
    if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
        return

    read_name = read.query_name
    if read_name not in reads_with_sa:
        reads_with_sa[read_name] = {}

    readp = reads_with_sa[read_name]
    if read.is_read1 and not read.is_supplementary:
        readp.setdefault("read1", []).append(read)
    elif read.is_read1 and read.is_supplementary:
        readp.setdefault("read1_s", []).append(read)
    elif read.is_read2 and not read.is_supplementary:
        readp.setdefault("read2", []).append(read)
    elif read.is_read2 and read.is_supplementary:
        readp.setdefault("read2_s", []).append(read)
    else:
        readp.setdefault("single", []).append(read)


def collect_ask_sa_reads_for_known_breakpoints(
    breakpoints: pd.DataFrame,
    bam_path: str | Path,
    window: int = 200,
    mapq: int = 0,
    nmmax: int = 3,
) -> dict:
    """
    Collect only reads near known breakpoint pairs, reducing ASK search space.
    """
    if pysam is None:
        raise ImportError("pysam is required for BAM-based ASK-search.")
    reads_with_sa: dict = {}
    bam_path = Path(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for bp in breakpoints.itertuples():
            for read in fetch_reads_around_breakpoint(
                bam,
                bp.Chrom1,
                bp.Coord1,
                bp.Chrom2,
                bp.Coord2,
                window,
            ):
                add_read_to_ask_sa_group(read, reads_with_sa, mapq=mapq, nmmax=nmmax)
    return reads_with_sa


def match_bp_duo_to_known_breakpoints(
    bp_duo: pd.DataFrame,
    breakpoints: pd.DataFrame,
    window: int = 200,
) -> pd.DataFrame:
    """
    Match ASK breakpoint-pair candidates to known ecDNA breakpoint pairs.
    """
    if bp_duo.empty:
        return pd.DataFrame()
    bps = infer_clip_from_known_breakpoints(breakpoints)
    rows: list[dict[str, object]] = []
    for _, cand in bp_duo.iterrows():
        obs1 = (cand["Chrom1"], cand["Coord1"], cand["Clip1"])
        obs2 = (cand["Chrom2"], cand["Coord2"], cand["Clip2"])
        for _, bp in bps.iterrows():
            if known_breakpoint_match(obs1, obs2, bp, window):
                rows.append(
                    {
                        "CircleID": bp["CircleID"],
                        "AmpliconID": bp.get("AmpliconID", ""),
                        "Filename": bp.get("Filename", ""),
                        "SampleID": bp.get("SampleID", ""),
                        "BreakpointID": bp["BreakpointID"],
                        "KnownChrom1": bp["Chrom1"],
                        "KnownCoord1": bp["Coord1"],
                        "KnownClip1": bp.get("Clip1", "R"),
                        "KnownChrom2": bp["Chrom2"],
                        "KnownCoord2": bp["Coord2"],
                        "KnownClip2": bp.get("Clip2", "L"),
                        "Gene": bp.get("Gene", ""),
                        "CancerGene": bp.get("CancerGene", ""),
                        "Chrom1": cand["Chrom1"],
                        "Coord1": cand["Coord1"],
                        "Clip1": cand["Clip1"],
                        "Chrom2": cand["Chrom2"],
                        "Coord2": cand["Coord2"],
                        "Clip2": cand["Clip2"],
                        "Count": cand.get("Count", 0),
                        "offset": cand.get("offset", None),
                        "Seq": cand.get("Seq", None),
                        "Readsname": cand.get("Readsname", []),
                    }
                )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).drop_duplicates(["CircleID", "BreakpointID", "Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2"])


def levenshtein_distance(a: object, b: object) -> int:
    """
    Compute Levenshtein edit distance for short junction sequences.
    """
    s1 = "" if pd.isna(a) else str(a)
    s2 = "" if pd.isna(b) else str(b)
    if s1 == s2:
        return 0
    if len(s1) < len(s2):
        s1, s2 = s2, s1
    previous = list(range(len(s2) + 1))
    for i, char1 in enumerate(s1, start=1):
        current = [i]
        for j, char2 in enumerate(s2, start=1):
            current.append(
                min(
                    previous[j] + 1,
                    current[j - 1] + 1,
                    previous[j - 1] + int(char1 != char2),
                )
            )
        previous = current
    return previous[-1]


def normalized_levenshtein(a: object, b: object) -> float:
    s1 = "" if pd.isna(a) else str(a)
    s2 = "" if pd.isna(b) else str(b)
    if not s1 and not s2:
        return 0.0
    return levenshtein_distance(s1.upper(), s2.upper()) / max(len(s1), len(s2), 1)


def junction_sequence_match(query_seq: object, ref_seq: object, max_norm_dist: float = 0.1) -> bool:
    if pd.isna(query_seq) or pd.isna(ref_seq):
        return False
    query = str(query_seq)
    ref = str(ref_seq)
    if not query or not ref or query in {"Known_Input", "PE_Support"} or ref in {"Known_Input", "PE_Support"}:
        return False
    return normalized_levenshtein(query, ref) <= max_norm_dist


def compute_jcs_tables(
    reference_junctions: pd.DataFrame,
    candidate_junctions: pd.DataFrame,
    window: int = 200,
    min_support_reads: int = 5,
    jcs_threshold: float = 0.5,
    max_norm_dist: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute Junction Concordance Score against the input reference structure.

    A reference junction is validated when a matched candidate junction has at
    least min_support_reads. If both reference and candidate junction sequences
    are available, normalized Levenshtein distance can additionally validate the
    sequence match.
    """
    if reference_junctions.empty:
        return pd.DataFrame(), pd.DataFrame()
    refs = infer_clip_from_known_breakpoints(reference_junctions)
    candidates = candidate_junctions.copy()
    rows: list[dict[str, object]] = []
    for _, ref in refs.iterrows():
        matches: list[dict[str, object]] = []
        for _, cand in candidates.iterrows():
            obs1 = (cand.get("Chrom1"), cand.get("Coord1"), cand.get("Clip1"))
            obs2 = (cand.get("Chrom2"), cand.get("Coord2"), cand.get("Clip2"))
            if not known_breakpoint_match(obs1, obs2, ref, window):
                continue
            count_value = cand.get("Count", 0)
            count = 0 if pd.isna(count_value) else int(float(count_value))
            seq_match = junction_sequence_match(cand.get("Seq", None), ref.get("RefJunctionSeq", None), max_norm_dist)
            matches.append(
                {
                    "Count": count,
                    "Seq": cand.get("Seq", ""),
                    "Readsname": cand.get("Readsname", []),
                    "SequenceMatch": seq_match,
                    "Chrom1": cand.get("Chrom1", ""),
                    "Coord1": cand.get("Coord1", ""),
                    "Clip1": cand.get("Clip1", ""),
                    "Chrom2": cand.get("Chrom2", ""),
                    "Coord2": cand.get("Coord2", ""),
                    "Clip2": cand.get("Clip2", ""),
                }
            )
        best = max(matches, key=lambda item: item["Count"], default={})
        support_count = int(best.get("Count", 0) or 0)
        validated = support_count >= min_support_reads
        rows.append(
            {
                "CircleID": ref.get("CircleID", ""),
                "AmpliconID": ref.get("AmpliconID", ""),
                "Filename": ref.get("Filename", ""),
                "SampleID": ref.get("SampleID", ""),
                "BreakpointID": ref.get("BreakpointID", ""),
                "Chrom1": ref.get("Chrom1", ""),
                "Coord1": ref.get("Coord1", ""),
                "Clip1": ref.get("Clip1", "R"),
                "Chrom2": ref.get("Chrom2", ""),
                "Coord2": ref.get("Coord2", ""),
                "Clip2": ref.get("Clip2", "L"),
                "SupportReadCount": support_count,
                "Validated": validated,
                "BestSeq": best.get("Seq", ""),
                "BestReadsname": best.get("Readsname", []),
                "SequenceMatch": best.get("SequenceMatch", False),
                "MatchedChrom1": best.get("Chrom1", ""),
                "MatchedCoord1": best.get("Coord1", ""),
                "MatchedClip1": best.get("Clip1", ""),
                "MatchedChrom2": best.get("Chrom2", ""),
                "MatchedCoord2": best.get("Coord2", ""),
                "MatchedClip2": best.get("Clip2", ""),
            }
        )
    junction_table = pd.DataFrame(rows)
    summary = (
        junction_table.groupby("CircleID", dropna=False, sort=False)
        .agg(
            total_reference_junctions=("BreakpointID", "nunique"),
            validated_junctions=("Validated", "sum"),
            total_support_reads=("SupportReadCount", "sum"),
        )
        .reset_index()
    )
    summary["JCS"] = summary["validated_junctions"] / summary["total_reference_junctions"]
    summary["Detected"] = summary["JCS"] > jcs_threshold
    summary["MinSupportReads"] = min_support_reads
    summary["JCSThreshold"] = jcs_threshold
    summary["MaxNormalizedLevenshtein"] = max_norm_dist
    return junction_table, summary


def search_known_bppairs_from_bam(
    breakpoints: pd.DataFrame,
    bam_path: str | Path,
    window: int = 200,
    mapq: int = 20,
    nmmax: int = 3,
    ask_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Targeted ASK-style breakpoint-pair search.

    The known ecDNA breakpoints restrict the reads scanned from BAM. The actual
    breakpoint-pair construction is still performed by ASK
    getbppair.get_bp_pair, so the returned table includes Count, offset, Seq,
    and Readsname.
    """
    getbppair = load_ask_getbppair(ask_path)
    bps = infer_clip_from_known_breakpoints(breakpoints)
    reads_with_sa = collect_ask_sa_reads_for_known_breakpoints(
        bps,
        bam_path=bam_path,
        window=window,
        mapq=0,
        nmmax=nmmax,
    )
    if reads_with_sa:
        bp_duo = getbppair.get_bp_pair(reads_with_sa, mapq=mapq)
    else:
        bp_duo = pd.DataFrame()
    matched = match_bp_duo_to_known_breakpoints(bp_duo, bps, window=window)
    return matched, bp_duo


def cigar_ref_length(cigar: str) -> int:
    """Reference-consuming length for a CIGAR string."""
    length = 0
    for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            length += int(size)
    return length


def sa_clip_side(sa: dict[str, object]) -> tuple[str, int] | tuple[None, None]:
    cigar = str(sa["cigar"])
    if not cigar:
        return None, None
    first = re.match(r"^(\d+)([SH])", cigar)
    last = re.search(r"(\d+)([SH])$", cigar)
    start = int(sa["pos"])
    end = start + cigar_ref_length(cigar) - 1
    if first:
        return "L", start
    if last:
        return "R", end
    return None, None


def coord_close(a: object, b: object, window: int) -> bool:
    return abs(int(a) - int(b)) <= window


def breakpoint_end_matches(
    observed_chrom: object,
    observed_coord: object,
    observed_clip: object,
    target_chrom: object,
    target_coord: object,
    target_clip: object,
    window: int,
) -> bool:
    return (
        str(observed_chrom) == str(target_chrom)
        and str(observed_clip) == str(target_clip)
        and coord_close(observed_coord, target_coord, window)
    )


def known_breakpoint_match(
    obs1: tuple[object, object, object],
    obs2: tuple[object, object, object],
    bp: pd.Series,
    window: int,
) -> bool:
    direct = breakpoint_end_matches(*obs1, bp["Chrom1"], bp["Coord1"], bp["Clip1"], window) and breakpoint_end_matches(
        *obs2, bp["Chrom2"], bp["Coord2"], bp["Clip2"], window
    )
    reverse = breakpoint_end_matches(*obs1, bp["Chrom2"], bp["Coord2"], bp["Clip2"], window) and breakpoint_end_matches(
        *obs2, bp["Chrom1"], bp["Coord1"], bp["Clip1"], window
    )
    return bool(direct or reverse)


def fetch_reads_around_breakpoint(bam, chrom1, coord1, chrom2, coord2, window: int):
    seen = set()
    for chrom, coord in [(chrom1, coord1), (chrom2, coord2)]:
        start = max(0, int(coord) - window)
        end = int(coord) + window + 1
        for read in bam.fetch(str(chrom), start, end):
            key = (read.query_name, read.flag, read.reference_id, read.reference_start)
            if key in seen:
                continue
            seen.add(key)
            yield read


def search_known_breakpoints_in_bam(
    breakpoints: pd.DataFrame,
    bam_path: str | Path,
    window: int = 200,
    mapq: int = 20,
    nmmax: int = 3,
    barcode_tags: Iterable[str] = ("CB", "CR", "BX", "BC"),
    include_normal: bool = True,
) -> pd.DataFrame:
    """
    Search a query BAM for evidence supporting known ecDNA breakpoint pairs.

    This is the targeted ASK-search mode. It follows ASK's evidence model:
    soft/hard-clipped primary reads and SA supplementary alignments define
    junction-supporting split reads, while normal reads spanning either genomic
    breakpoint side can be collected as chromosomal/reference evidence.
    """
    if pysam is None:
        raise ImportError("pysam is required for BAM-based ASK-search.")

    bps = infer_clip_from_known_breakpoints(breakpoints)
    require_columns(bps, ["CircleID", "BreakpointID", "Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2"])

    rows: list[dict[str, object]] = []
    bam_path = Path(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for _, bp in bps.iterrows():
            for read in fetch_reads_around_breakpoint(bam, bp["Chrom1"], bp["Coord1"], bp["Chrom2"], bp["Coord2"], window):
                if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
                    continue

                clip_side, clip_coord = read_clip_side(read)
                primary_obs = None
                if clip_side is not None:
                    primary_obs = (read.reference_name, clip_coord, clip_side)

                matched = False
                for sa in parse_sa_tag(read):
                    if sa["mapq"] < mapq or sa["nm"] > nmmax:
                        continue
                    sa_side, sa_coord = sa_clip_side(sa)
                    if primary_obs is None or sa_side is None:
                        continue
                    sa_obs = (sa["chrom"], sa_coord, sa_side)
                    if known_breakpoint_match(primary_obs, sa_obs, bp, window):
                        rows.append(
                            bam_evidence_row(
                                bp,
                                bam_path,
                                read,
                                "split_SA",
                                read_barcode(read, barcode_tags),
                                mate_chrom=sa["chrom"],
                                mate_coord=sa_coord,
                                mate_clip=sa_side,
                                mate_cigar=sa["cigar"],
                                distance_to_bp1=min(abs(int(primary_obs[1]) - int(bp["Coord1"])), abs(int(sa_obs[1]) - int(bp["Coord1"]))),
                                distance_to_bp2=min(abs(int(primary_obs[1]) - int(bp["Coord2"])), abs(int(sa_obs[1]) - int(bp["Coord2"]))),
                            )
                        )
                        matched = True
                        break

                if matched:
                    continue

                if primary_obs is not None:
                    side1 = breakpoint_end_matches(*primary_obs, bp["Chrom1"], bp["Coord1"], bp["Clip1"], window)
                    side2 = breakpoint_end_matches(*primary_obs, bp["Chrom2"], bp["Coord2"], bp["Clip2"], window)
                    if side1 or side2:
                        rows.append(
                            bam_evidence_row(
                                bp,
                                bam_path,
                                read,
                                "single_clip",
                                read_barcode(read, barcode_tags),
                                distance_to_bp1=abs(int(primary_obs[1]) - int(bp["Coord1"])),
                                distance_to_bp2=abs(int(primary_obs[1]) - int(bp["Coord2"])),
                            )
                        )
                        continue

                if include_normal:
                    spans_bp1 = str(read.reference_name) == str(bp["Chrom1"]) and read.reference_start <= int(bp["Coord1"]) <= read.reference_end
                    spans_bp2 = str(read.reference_name) == str(bp["Chrom2"]) and read.reference_start <= int(bp["Coord2"]) <= read.reference_end
                    if spans_bp1 or spans_bp2:
                        rows.append(
                            bam_evidence_row(
                                bp,
                                bam_path,
                                read,
                                "normal_spanning",
                                read_barcode(read, barcode_tags),
                                breakpoint_side="Bp1" if spans_bp1 else "Bp2",
                            )
                        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "CircleID",
                "AmpliconID",
                "Filename",
                "SampleID",
                "BreakpointID",
                "Chrom1",
                "Coord1",
                "Clip1",
                "Chrom2",
                "Coord2",
                "Clip2",
                "Gene",
                "CancerGene",
                "Bam",
                "Readname",
                "Barcode",
                "support_type",
                "MAPQ",
                "NM",
                "CIGAR",
            ]
        )
    return pd.DataFrame(rows).drop_duplicates()


def bam_evidence_row(
    bp: pd.Series,
    bam_path: Path,
    read,
    support_type: str,
    barcode: object,
    mate_chrom: object = None,
    mate_coord: object = None,
    mate_clip: object = None,
    mate_cigar: object = None,
    distance_to_bp1: object = None,
    distance_to_bp2: object = None,
    breakpoint_side: object = None,
) -> dict[str, object]:
    clip_side, clip_coord = read_clip_side(read)
    return {
        "CircleID": bp["CircleID"],
        "AmpliconID": bp.get("AmpliconID", ""),
        "Filename": bp.get("Filename", ""),
        "SampleID": bp.get("SampleID", ""),
        "BreakpointID": bp["BreakpointID"],
        "Chrom1": bp["Chrom1"],
        "Coord1": bp["Coord1"],
        "Clip1": bp.get("Clip1", "R"),
        "Chrom2": bp["Chrom2"],
        "Coord2": bp["Coord2"],
        "Clip2": bp.get("Clip2", "L"),
        "Gene": bp.get("Gene", ""),
        "CancerGene": bp.get("CancerGene", ""),
        "Bam": str(bam_path),
        "Readname": read.query_name,
        "Barcode": barcode,
        "support_type": support_type,
        "breakpoint_side": breakpoint_side,
        "ReadChrom": read.reference_name,
        "ReadCoord": clip_coord if clip_coord is not None else read.reference_start,
        "ReadClip": clip_side,
        "MateChrom": mate_chrom,
        "MateCoord": mate_coord,
        "MateClip": mate_clip,
        "MAPQ": read.mapping_quality,
        "NM": safe_get_tag(read, "NM"),
        "CIGAR": read.cigarstring,
        "MateCIGAR": mate_cigar,
        "is_reverse": read.is_reverse,
        "is_supplementary": read.is_supplementary,
        "distance_to_bp1": distance_to_bp1,
        "distance_to_bp2": distance_to_bp2,
    }


def evidence_from_targeted_bppair_matches(
    matched_pairs: pd.DataFrame,
    bam_path: str | Path,
    barcode_lookup: dict[str, object] | None = None,
) -> pd.DataFrame:
    """
    Expand matched ASK bp-pair rows into read-level evidence.
    """
    rows: list[dict[str, object]] = []
    if matched_pairs.empty:
        return pd.DataFrame()
    for _, row in matched_pairs.iterrows():
        readnames = list(dict.fromkeys(parse_list_cell(row.get("Readsname", []))))
        if not readnames:
            readnames = [None]
        for readname in readnames:
            rows.append(
                {
                    "CircleID": row["CircleID"],
                    "AmpliconID": row.get("AmpliconID", ""),
                    "Filename": row.get("Filename", ""),
                    "SampleID": row.get("SampleID", ""),
                    "BreakpointID": row["BreakpointID"],
                    "Chrom1": row["KnownChrom1"],
                    "Coord1": row["KnownCoord1"],
                    "Clip1": row["KnownClip1"],
                    "Chrom2": row["KnownChrom2"],
                    "Coord2": row["KnownCoord2"],
                    "Clip2": row["KnownClip2"],
                    "Gene": row.get("Gene", ""),
                    "CancerGene": row.get("CancerGene", ""),
                    "Bam": str(bam_path),
                    "Readname": readname,
                    "Barcode": barcode_lookup.get(readname) if barcode_lookup else None,
                    "support_type": "split_SA",
                    "CandidateChrom1": row["Chrom1"],
                    "CandidateCoord1": row["Coord1"],
                    "CandidateClip1": row["Clip1"],
                    "CandidateChrom2": row["Chrom2"],
                    "CandidateCoord2": row["Coord2"],
                    "CandidateClip2": row["Clip2"],
                    "Count": row.get("Count", None),
                    "offset": row.get("offset", None),
                    "Seq": row.get("Seq", None),
                }
            )
    return pd.DataFrame(rows)


def collect_read_barcodes_from_bam(
    bam_path: str | Path,
    readnames: Iterable[str],
    barcode_tags: Iterable[str] = ("CB", "CR", "BX", "BC"),
) -> dict[str, object]:
    """
    Resolve read names to barcode tags by scanning the BAM once.
    """
    if pysam is None:
        raise ImportError("pysam is required for BAM-based ASK-search.")
    wanted = {str(readname) for readname in readnames if readname is not None and not pd.isna(readname)}
    if not wanted:
        return {}
    out = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.query_name in wanted and read.query_name not in out:
                out[read.query_name] = read_barcode(read, barcode_tags)
                if len(out) == len(wanted):
                    break
    return out


def genome_bin_count_from_bam(
    bam_path: str | Path,
    binsize: int = 10000,
    mapq: int = 20,
    nmmax: int = 3,
) -> pd.DataFrame:
    """
    ASK-style genome bin read count and CN normalization.

    CN is estimated as 2 * bin_count / mean(bin_count), matching the simple
    normalization used by ASK covbam.region_count.
    """
    if pysam is None:
        raise ImportError("pysam is required for BAM-based copy-number estimation.")
    counts: dict[tuple[str, int], int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
                continue
            coord = int(read.reference_start // binsize * binsize)
            key = (read.reference_name, coord)
            counts[key] = counts.get(key, 0) + 1
    if not counts:
        return pd.DataFrame(columns=["Chrom", "Coord", "Count", "CN"])
    df = pd.DataFrame(
        [(chrom, coord, count) for (chrom, coord), count in counts.items()],
        columns=["Chrom", "Coord", "Count"],
    ).sort_values(["Chrom", "Coord"])
    mean_count = df["Count"].mean()
    df["CN"] = 2 * df["Count"] / mean_count if mean_count else 0
    return df.reset_index(drop=True)


def count_reads_in_interval(
    bam,
    chrom: object,
    start: object,
    end: object,
    mapq: int = 20,
    nmmax: int = 3,
) -> int:
    """
    Count filtered reads overlapping an interval.

    ASK circular tables are 1-based-like genomic coordinates; pysam fetch is
    0-based half-open, so start is shifted down by one.
    """
    start0 = max(0, int(start) - 1)
    end0 = int(end)
    seen = set()
    for read in bam.fetch(str(chrom), start0, end0):
        if not ask_check_read(read, mapq=mapq, nmmax=nmmax):
            continue
        key = (read.query_name, read.flag, read.reference_id, read.reference_start)
        seen.add(key)
    return len(seen)


def weighted_cn_for_interval(bin_count: pd.DataFrame, chrom: object, start: object, end: object, binsize: int) -> tuple[float | None, int]:
    """
    Weighted mean CN over bins overlapping one segment.
    """
    if bin_count.empty:
        return None, 0
    start0 = max(0, int(start) - 1)
    end0 = int(end)
    sub = bin_count[
        (bin_count["Chrom"].astype(str) == str(chrom))
        & (bin_count["Coord"] < end0)
        & ((bin_count["Coord"] + binsize) > start0)
    ].copy()
    if sub.empty:
        return None, 0
    overlaps = []
    for row in sub.itertuples():
        bin_start = int(row.Coord)
        bin_end = bin_start + binsize
        overlaps.append(max(0, min(end0, bin_end) - max(start0, bin_start)))
    sub["OverlapBp"] = overlaps
    total = sub["OverlapBp"].sum()
    if total <= 0:
        return None, int(sub.shape[0])
    cn = float((sub["CN"] * sub["OverlapBp"]).sum() / total)
    return cn, int(sub.shape[0])


def estimate_segment_copy_number(
    circular_df: pd.DataFrame,
    bam_path: str | Path,
    binsize: int = 10000,
    mapq: int = 20,
    nmmax: int = 3,
    bin_count: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Estimate copy number for each known ecDNA segment in the query BAM.

    Returns segment-level CN table and the genome bin-count table used for
    normalization.
    """
    require_columns(circular_df, ["CircleID", "AmpliconID", "Filename", "Chrom", "Start", "End", "Strand"])
    if bin_count is None:
        bin_count = genome_bin_count_from_bam(bam_path, binsize=binsize, mapq=mapq, nmmax=nmmax)

    rows: list[dict[str, object]] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for idx, seg in circular_df.reset_index(drop=True).iterrows():
            direct_count = count_reads_in_interval(
                bam,
                seg["Chrom"],
                seg["Start"],
                seg["End"],
                mapq=mapq,
                nmmax=nmmax,
            )
            length = max(1, int(seg["End"]) - int(seg["Start"]) + 1)
            cn, n_bins = weighted_cn_for_interval(bin_count, seg["Chrom"], seg["Start"], seg["End"], binsize)
            rows.append(
                {
                    "CircleID": seg["CircleID"],
                    "AmpliconID": seg["AmpliconID"],
                    "Filename": seg["Filename"],
                    "SegmentIndex": idx,
                    "Chrom": seg["Chrom"],
                    "Start": int(seg["Start"]),
                    "End": int(seg["End"]),
                    "Strand": seg["Strand"],
                    "Length": length,
                    "Gene": seg.get("Gene", ""),
                    "CancerGene": seg.get("CancerGene", ""),
                    "ReadCount": direct_count,
                    "ReadsPerKb": direct_count / (length / 1000),
                    "CN": cn,
                    "OverlappingBins": n_bins,
                    "Binsize": binsize,
                }
            )
    return pd.DataFrame(rows), bin_count


def summarize_circle_copy_number(segment_cn: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize segment copy-number estimates at circle level.
    """
    if segment_cn.empty:
        return pd.DataFrame()
    df = segment_cn.copy()
    summary = (
        df.groupby(["CircleID", "AmpliconID", "Filename"], dropna=False)
        .agg(
            Chrom=("Chrom", lambda x: ";".join(sorted(set(map(str, x))))),
            Start=("Start", "min"),
            End=("End", "max"),
            Seg_num=("SegmentIndex", "count"),
            Length=("Length", "sum"),
            ReadCount_sum=("ReadCount", "sum"),
            ReadCount_mean=("ReadCount", "mean"),
            ReadsPerKb_mean=("ReadsPerKb", "mean"),
            CN_mean=("CN", "mean"),
            CN_median=("CN", "median"),
            CN_std=("CN", "std"),
            CN_min=("CN", "min"),
            CN_max=("CN", "max"),
            Gene_num=("Gene", lambda x: len(set(g for v in x.dropna() for g in split_genes(v)))),
            Cancergene_num=("CancerGene", lambda x: len(set(g for v in x.dropna() for g in split_genes(v)))),
        )
        .reset_index()
    )
    summary["CN_std"] = summary["CN_std"].fillna(0)
    return summary


def load_read_barcode_tables(paths: Iterable[str | Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        path = Path(path)
        if path.is_dir():
            files = sorted(list(path.glob("*.csv")) + list(path.glob("*.tsv")))
        else:
            files = [path]
        for file in files:
            df = read_table(file)
            if "Filename" not in df.columns:
                df["Filename"] = file.name
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["ReadID", "Barcode", "Filename"])
    out = pd.concat(frames, ignore_index=True)
    require_columns(out, ["ReadID", "Barcode"])
    return out


def find_junction_evidence(
    breakpoints: pd.DataFrame,
    bp_pair: pd.DataFrame,
    read_barcodes: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Match circular breakpoint pairs to ASK breakpoint-pair rows and expand
    read-level support into a tidy evidence table.
    """
    require_columns(breakpoints, ["BreakpointID", "CircleID", "Chrom1", "Coord1", "Chrom2", "Coord2"])
    require_columns(bp_pair, ["Chrom1", "Coord1", "Chrom2", "Coord2", "Readname"])
    bp_pair = add_breakpoint_key(bp_pair)
    breakpoints = add_breakpoint_key(breakpoints)

    read_lookup = None
    if read_barcodes is not None and not read_barcodes.empty:
        read_lookup = read_barcodes.groupby("ReadID", sort=False)

    rows: list[dict[str, object]] = []
    bp_groups = {key: group for key, group in bp_pair.groupby("_bp_key", sort=False)}
    for _, bp in breakpoints.iterrows():
        matches = bp_groups.get(bp["_bp_key"])
        if matches is None:
            continue

        for match in matches.itertuples():
            readnames = parse_list_cell(getattr(match, "Readname", None))
            readbarcodes = parse_list_cell(getattr(match, "Readbarcode", None))
            if len(readbarcodes) < len(readnames):
                readbarcodes.extend([None] * (len(readnames) - len(readbarcodes)))

            for read_idx, readname in enumerate(readnames):
                ask_barcode = readbarcodes[read_idx] if read_idx < len(readbarcodes) else None
                barcode_rows = None
                if read_lookup is not None and readname in read_lookup.groups:
                    barcode_rows = read_lookup.get_group(readname)

                if barcode_rows is None or barcode_rows.empty:
                    rows.append(evidence_row(bp, readname, ask_barcode, None, None))
                    continue

                for barcode_row in barcode_rows.itertuples():
                    rows.append(
                        evidence_row(
                            bp,
                            readname,
                            ask_barcode,
                            getattr(barcode_row, "Barcode", None),
                            getattr(barcode_row, "Filename", None),
                        )
                    )

    return pd.DataFrame(rows)


def evidence_row(bp: pd.Series, readname: object, ask_barcode: object, cell_barcode: object, read_filename: object) -> dict[str, object]:
    return {
        "CircleID": bp["CircleID"],
        "AmpliconID": bp["AmpliconID"],
        "Filename": bp["Filename"],
        "SampleID": bp["SampleID"],
        "BreakpointID": bp["BreakpointID"],
        "Chrom1": bp["Chrom1"],
        "Coord1": bp["Coord1"],
        "Chrom2": bp["Chrom2"],
        "Coord2": bp["Coord2"],
        "Gene": bp.get("Gene", ""),
        "CancerGene": bp.get("CancerGene", ""),
        "Readname": readname,
        "ASKReadBarcode": ask_barcode,
        "Barcode": cell_barcode,
        "ReadFilename": read_filename,
    }


def normalize_gene_label(gene_text: object, target_genes: Iterable[str] = TARGET_GENES) -> str:
    genes = set()
    text = "" if pd.isna(gene_text) else str(gene_text)
    for gene in target_genes:
        if re.search(rf"(^|[;\s,_]){re.escape(gene)}($|[;\s,_])", text):
            genes.add(gene)
    return "_".join(sorted(genes)) if genes else "ecDNA"


def make_barcode_id(row, sample_prefix: str | None = None) -> str:
    barcode = getattr(row, "Barcode", None)
    if pd.isna(barcode) or barcode is None:
        return ""
    if sample_prefix:
        sample = getattr(row, "ReadFilename", None) or getattr(row, "SampleID", None)
        return f"{sample_prefix}_{sample}_{barcode}"
    return str(barcode)


def build_matrices(
    evidence: pd.DataFrame,
    min_reads: int = 1,
    min_breakpoints: int = 1,
    sample_prefix: str | None = None,
    support_types: Iterable[str] = ("split_SA", "single_clip"),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Return:
      cell_breakpoint_matrix: cells x breakpoint IDs, values are read counts.
      cell_circle_summary: one row per cell/circle.
      cell_annotation: one row per cell with ecDNA+/No classification.
    """
    if evidence.empty:
        empty = pd.DataFrame()
        return empty, empty, empty

    df = evidence.copy()
    if "support_type" in df.columns and support_types is not None:
        df = df[df["support_type"].isin(set(support_types))].copy()
    df = df[df["Barcode"].notna()].copy()
    if df.empty:
        empty = pd.DataFrame()
        return empty, empty, empty

    df["CellID"] = [make_barcode_id(row, sample_prefix) for row in df.itertuples()]
    df = df[df["CellID"] != ""]
    df["GeneLabel"] = df["CancerGene"].map(normalize_gene_label)

    dedup = df.drop_duplicates(["CellID", "CircleID", "BreakpointID", "Readname"])
    cell_breakpoint_matrix = (
        dedup.groupby(["CellID", "BreakpointID"], sort=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    cell_circle_summary = (
        dedup.groupby(["CellID", "CircleID", "GeneLabel"], sort=False)
        .agg(
            n_reads=("Readname", "nunique"),
            n_breakpoints=("BreakpointID", "nunique"),
            first_sample=("SampleID", "first"),
            first_filename=("Filename", "first"),
        )
        .reset_index()
    )
    cell_circle_summary["is_ecDNA_positive"] = (
        (cell_circle_summary["n_reads"] >= min_reads)
        & (cell_circle_summary["n_breakpoints"] >= min_breakpoints)
    )

    positive = cell_circle_summary[cell_circle_summary["is_ecDNA_positive"]].copy()
    if positive.empty:
        cell_annotation = pd.DataFrame({"CellID": sorted(df["CellID"].unique()), "ecDNA_status": "No"})
    else:
        cell_annotation = (
            positive.groupby("CellID", sort=False)
            .agg(
                ecDNA_status=("GeneLabel", lambda x: "_".join(sorted(set(x)))),
                n_ecDNA_circles=("CircleID", "nunique"),
                n_ecDNA_breakpoints=("n_breakpoints", "sum"),
                n_ecDNA_reads=("n_reads", "sum"),
            )
            .reset_index()
        )
        all_cells = pd.DataFrame({"CellID": sorted(df["CellID"].unique())})
        cell_annotation = all_cells.merge(cell_annotation, on="CellID", how="left")
        cell_annotation["ecDNA_status"] = cell_annotation["ecDNA_status"].fillna("No")
        for col in ["n_ecDNA_circles", "n_ecDNA_breakpoints", "n_ecDNA_reads"]:
            cell_annotation[col] = cell_annotation[col].fillna(0).astype(int)

    return cell_breakpoint_matrix, cell_circle_summary, cell_annotation


def run_extract(args: argparse.Namespace) -> None:
    samples = [parse_sample_arg(value) for value in args.circular]
    circular = load_circular_tables(samples)
    if args.target_genes:
        target_genes = set(args.target_genes.split(","))
        keep = circular.groupby("CircleID")["CancerGene"].transform(
            lambda values: any(gene in target_genes for value in values.dropna() for gene in split_genes(value))
        )
        circular = circular[keep].copy()
    breakpoints = extract_breakpoint_pairs(circular)
    write_table(breakpoints, args.output)


def run_evidence(args: argparse.Namespace) -> None:
    breakpoints = read_table(args.breakpoints)
    bp_pair = read_table(args.bp_pair)
    read_barcodes = load_read_barcode_tables(args.read_barcodes or [])
    evidence = find_junction_evidence(breakpoints, bp_pair, read_barcodes)
    write_table(evidence, args.output)


def run_classify(args: argparse.Namespace) -> None:
    evidence = read_table(args.evidence)
    matrix, summary, annotation = build_matrices(
        evidence,
        min_reads=args.min_reads,
        min_breakpoints=args.min_breakpoints,
        sample_prefix=args.sample_prefix,
        support_types=args.support_types.split(",") if args.support_types else None,
    )
    outdir = Path(args.outdir)
    write_table(matrix, outdir / "cell_by_breakpoint_matrix.tsv")
    write_table(summary, outdir / "cell_circle_summary.tsv")
    write_table(annotation, outdir / "cell_ecDNA_annotation.tsv")


def run_all(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = [parse_sample_arg(value) for value in args.circular]
    circular = load_circular_tables(samples)
    if args.target_genes:
        target_genes = set(args.target_genes.split(","))
        keep = circular.groupby("CircleID")["CancerGene"].transform(
            lambda values: any(gene in target_genes for value in values.dropna() for gene in split_genes(value))
        )
        circular = circular[keep].copy()

    breakpoints = extract_breakpoint_pairs(circular)
    write_table(breakpoints, outdir / "ecDNA_breakpoint_pairs.tsv")

    bp_pair = read_table(args.bp_pair)
    read_barcodes = load_read_barcode_tables(args.read_barcodes or [])
    evidence = find_junction_evidence(breakpoints, bp_pair, read_barcodes)
    write_table(evidence, outdir / "junction_read_evidence.tsv")

    matrix, summary, annotation = build_matrices(
        evidence,
        min_reads=args.min_reads,
        min_breakpoints=args.min_breakpoints,
        sample_prefix=args.sample_prefix,
    )
    write_table(matrix, outdir / "cell_by_breakpoint_matrix.tsv")
    write_table(summary, outdir / "cell_circle_summary.tsv")
    write_table(annotation, outdir / "cell_ecDNA_annotation.tsv")


def run_bam_search(args: argparse.Namespace) -> None:
    if pysam is None:
        raise ImportError("pysam is required for BAM-based ASK-search.")

    samples = [parse_sample_arg(value) for value in args.circular]
    circular = load_circular_tables(samples)
    if args.target_genes:
        target_genes = set(args.target_genes.split(","))
        keep = circular.groupby("CircleID")["CancerGene"].transform(
            lambda values: any(gene in target_genes for value in values.dropna() for gene in split_genes(value))
        )
        circular = circular[keep].copy()

    breakpoints = extract_breakpoint_pairs(circular)
    breakpoints = infer_clip_from_known_breakpoints(breakpoints)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ask_outputs = ask_style_output_paths(outdir, args.bam, args.outprefix)
    cleanup_legacy_bam_search_outputs(outdir)
    genefile, sefile, cgfile = resolve_annotation_files(
        genome=args.genome,
        genefile=args.genefile,
        sefile=args.sefile,
        cgfile=args.cgfile,
    )
    write_table(circular, outdir / "known_ecDNA_segments.tsv")
    write_table(breakpoints, outdir / "known_ecDNA_breakpoint_pairs.tsv")

    target_chroms = target_chromosomes_from_circular(circular)
    bp_all, clip_bg, bin_count, bp_duo_bam, reads_with_sa = target_process_alignment(
        bam_path=args.bam,
        breakpoints=breakpoints,
        target_chroms=target_chroms,
        circular_df=circular,
        cn_background=args.cn_background,
        window=args.window,
        binsize=args.binsize,
        mapq=args.mapq,
        nmmax=args.nmmax,
        ask_path=args.ask_path,
    )
    write_table(bin_count, ask_outputs["output_bincnt"])
    clip_bg.to_csv(ask_outputs["output_clip"], sep="\t", index=False, header=False)
    write_pickle([bp_all, clip_bg, bin_count, bp_duo_bam, str(args.bam)], ask_outputs["output_pdat_1"])

    cn_amp, cn_seg, bin_norm = target_detect_amplified_segment(
        bin_count=bin_count,
        bp_all=bp_all,
        circular_df=circular,
        binsize=args.binsize,
        ask_path=args.ask_path,
        min_cn=args.min_cn,
        min_segsize_scale=args.min_segsize_scale,
    )
    cn_amp = annotate_amplicon_table(
        cn_amp,
        genefile=genefile,
        sefile=None,
        cgfile=cgfile,
        ask_path=args.ask_path,
    )
    write_table(cn_amp, ask_outputs["output_cnamp"])
    write_table(cn_seg, ask_outputs["output_cnseg"])
    write_table(bin_norm, ask_outputs["output_binnorm"])
    write_pickle([cn_amp, cn_seg, bin_norm], ask_outputs["output_pdat_2"])

    known_seed = known_breakpoint_seed_table(breakpoints)
    write_table(known_seed, outdir / "known_breakpoint_seed.tsv")
    bp_duo, bp_cand_stats, bp_duo_sa_near_known = target_detect_bp_pair(
        bam_path=args.bam,
        known_seed=known_seed,
        bp_duo_bam=bp_duo_bam,
        bp_all=bp_all,
        cn_amp=cn_amp,
        binsize=args.binsize,
        window=args.window,
        min_nt=args.min_nt,
        mapq=args.mapq,
        nmmax=args.nmmax,
        ask_path=args.ask_path,
        include_sa_candidates=not args.no_sa_candidates,
    )
    write_table(bp_cand_stats, ask_outputs["output_bpcand"])
    write_table(bp_duo, ask_outputs["output_bpduo"])
    write_table(bp_duo, ask_outputs["output_align"])
    _, jcs_summary = compute_jcs_tables(
        breakpoints,
        bp_duo,
        window=args.window,
        min_support_reads=args.jcs_min_support,
        jcs_threshold=args.min_jcs,
        max_norm_dist=args.jcs_max_norm_dist,
    )
    write_table(jcs_summary, ask_outputs["output_jcs"])
    write_pickle([bp_duo, bp_cand_stats], ask_outputs["output_pdat_3"])

    circ_anno, line_anno, bp_pair, seg, circ_score, bp_fine = target_construct_amplicon(
        bp_duo=bp_duo,
        bp_cand_stats=bp_cand_stats,
        cn_amp=cn_amp,
        bin_norm=bin_norm if not bin_norm.empty else bin_count,
        known_single_segments=known_single_segment_circles(circular),
        min_junc_cnt=args.min_junc_cnt,
        bpp_min_dist=args.bpp_min_dist,
        segment_restrict=args.segment_restrict,
        subseg=not args.no_subseg,
        rm_amp_df=args.rm_amp_df,
        ask_path=args.ask_path,
    )
    circ_anno = annotate_amplicon_table(
        circ_anno,
        genefile=genefile,
        sefile=sefile,
        cgfile=cgfile,
        ask_path=args.ask_path,
    )
    line_anno = annotate_amplicon_table(
        line_anno,
        genefile=genefile,
        sefile=sefile,
        cgfile=cgfile,
        ask_path=args.ask_path,
    )
    if not circ_anno.empty and not (bin_norm if not bin_norm.empty else bin_count).empty:
        try:
            ggraph = load_ask_modules(args.ask_path)["ggraph"]
            circ_score = ggraph.add_stats_circ(
                circ_anno,
                bin_norm if not bin_norm.empty else bin_count,
                binsize=args.binsize,
            )
        except Exception:
            pass
    write_table(bp_pair, ask_outputs["output_bppair"])
    write_table(seg, ask_outputs["output_seg"])
    write_table(circ_anno, ask_outputs["output_circ"])
    write_table(line_anno, ask_outputs["output_line"])
    write_table(circ_score, ask_outputs["output_circ_stat"])
    write_junction_sequence_outputs(
        bp_pair,
        bam_path=args.bam,
        output_align=ask_outputs["output_align"],
        output_align_dir=ask_outputs["output_align_dir"],
        ask_path=args.ask_path,
        circ_anno=circ_anno,
    )
    plot_target_amplicons(
        circ_anno,
        line_anno,
        cn_amp,
        genefile=genefile,
        sefile=sefile,
        bin_norm=bin_norm if not bin_norm.empty else bin_count,
        binsize=args.binsize,
        fig_dir=ask_outputs["outfig"],
        ask_path=args.ask_path,
    )
    write_pickle([circ_anno, line_anno, bp_pair, seg], ask_outputs["output_pdat_4"])
    write_ask_stats(ask_outputs["output_stats"], cn_amp, bp_cand_stats, bp_duo, circ_anno, line_anno)

    cleanup_legacy_bam_search_outputs(outdir)


def add_bam_search_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument("--circular", action="append", required=True, help="sample_id[,filename]=known_ecDNA_structure.tsv")
    parser.add_argument("--bam", required=True, help="Query BAM/CRAM indexed for random access.")
    parser.add_argument("--target-genes", default=None)
    parser.add_argument("--genome", default="hg38", help="Genome build for default ASK annotations, e.g. hg19, hg38, mm10.")
    parser.add_argument("--genefile", default=None, help="Optional BED12 gene annotation file.")
    parser.add_argument("--sefile", default=None, help="Optional super-enhancer BED file.")
    parser.add_argument("--cgfile", default=None, help="Optional cancer gene census file.")
    parser.add_argument("--window", type=int, default=200)
    parser.add_argument("--binsize", type=int, default=10000)
    parser.add_argument(
        "--cn-background",
        choices=["genome", "target"],
        default="genome",
        help="Background bins used for CN normalization. Use genome to match ASK-style plots.",
    )
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--nmmax", type=int, default=1)
    parser.add_argument("--min-cn", type=float, default=None, help="CN cutoff for targeted amplified segment calling.")
    parser.add_argument("--min-segsize-scale", type=int, default=1, help="Minimum amplified segment size in bins.")
    parser.add_argument("--min-nt", type=int, default=5, help="Minimum matched nucleotides for ASK bpjoint.get_bppair.")
    parser.add_argument("--min-junc-cnt", type=int, default=1, help="Minimum junction read count used before DFS construction.")
    parser.add_argument("--bpp-min-dist", type=int, default=50, help="Minimum same-chromosome distance for breakpoint-pair filtering.")
    parser.add_argument("--jcs-min-support", type=int, default=5, help="Minimum supporting reads required to validate one reference junction.")
    parser.add_argument("--min-jcs", type=float, default=0.5, help="Circle-level Junction Concordance Score detection threshold.")
    parser.add_argument("--jcs-max-norm-dist", type=float, default=0.1, help="Maximum normalized Levenshtein distance for optional junction-sequence matching.")
    parser.add_argument("--barcode-tags", default="CB,CR,BX,BC")
    parser.add_argument("--ask-path", dest="ask_path", default=None, help="Path to ASK source directory containing getbppair.py.")
    parser.add_argument("--no-normal", action="store_true", help="Do not collect normal reads spanning breakpoint positions.")
    parser.add_argument("--no-sa-candidates", action="store_true", help="Use only known breakpoint seeds, not nearby SA-derived candidates.")
    parser.add_argument("--segment-restrict", action="store_true", help="Restrict ASK segment inference to clean breakpoints.")
    parser.add_argument("--no-subseg", action="store_true", help="Disable ASK subsegment inference inside CN amplicons.")
    parser.add_argument("--rm-amp-df", action="store_true", help="Remove amplicon boundary breakpoints after ASK subsegment inference.")
    parser.add_argument("--min-reads", type=int, default=1)
    parser.add_argument("--min-breakpoints", type=int, default=2)
    parser.add_argument("--support-types", default="split_SA,single_clip")
    parser.add_argument("--sample-prefix", default=None)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--outprefix", default=None, help="Optional ASK-style output prefix. Defaults to outdir/<bam-stem>.")
    parser.set_defaults(func=run_bam_search)
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ASK-search: targeted search for known ecDNA breakpoint evidence directly from a query BAM/CRAM."
    )
    return add_bam_search_arguments(parser)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
