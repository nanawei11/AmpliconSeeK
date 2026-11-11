import ast
import json
from collections import Counter
from pathlib import Path

import pandas as pd

try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None


BARCODE_TAGS = ("CB", "CR", "BX", "BC")


def parse_list_cell(value):
    if isinstance(value, list):
        return value
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text in {"[]", "None", "nan"}:
        return []
    try:
        parsed = ast.literal_eval(text)
    except Exception:
        return [text]
    if isinstance(parsed, list):
        return parsed
    return [parsed]


def clean_barcode(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text in {"None", "nan"}:
        return None
    return text


def clean_barcode_list(value):
    return [barcode for barcode in (clean_barcode(x) for x in parse_list_cell(value)) if barcode]


def clean_barcode_columns(df):
    if df.empty:
        return df
    out = df.copy()
    for col in ("Readbarcode", "Readsbarcode", "BestReadsbarcode"):
        if col in out.columns:
            out[col] = out[col].map(clean_barcode_list)
    return out


def has_single_cell_barcodes(df):
    if df.empty:
        return False
    for col in ("Readbarcode", "Readsbarcode", "BestReadsbarcode"):
        if col in df.columns and df[col].map(lambda value: len(clean_barcode_list(value)) > 0).any():
            return True
    return False


def read_barcode(read, tags=BARCODE_TAGS):
    for tag in tags:
        if read.has_tag(tag):
            return clean_barcode(read.get_tag(tag))
    return None


def junction_id(chrom1, coord1, clip1, chrom2, coord2, clip2):
    return f"{chrom1}:{int(coord1)}:{clip1}|{chrom2}:{int(coord2)}:{clip2}"


def _segment_exit_breakend(row):
    """
    Segment 'exit' breakend (where the circle leaves this segment for the next).
    Clip is the breakpoint's intrinsic type: R=segment right end (End),
    L=segment left end (Start). The exit depends on strand:
      forward (+): (End, R)    reverse (-): (Start, L)
    """
    if str(row.get("Strand", "+")) == "-":
        return row["Chrom"], int(row["Start"]), "L"
    return row["Chrom"], int(row["End"]), "R"


def _segment_entry_breakend(row):
    """
    Segment 'entry' breakend (where the circle enters this segment from the
    previous one). Depends on strand:
      forward (+): (Start, L)    reverse (-): (End, R)
    """
    if str(row.get("Strand", "+")) == "-":
        return row["Chrom"], int(row["End"]), "R"
    return row["Chrom"], int(row["Start"]), "L"


def infer_circular_junctions(circular_df):
    """
    Infer junctions (breakpoint pairs) between adjacent segments of a circle.

    FIX(20260717): clip is taken from segment strand instead of hardcoding
      (End,R)->(Start,L). With inverted (-) segments the intrinsic junction clip
      is RR/LL, so hardcoded RL would not match bp_pair and would zero out that
      junction's support/normal (same root cause as in ecDNA_search).
      junction = current segment exit -> next segment entry; a single-segment
      circle is its own exit -> entry (closure).
    """
    rows = []
    required = {"AmpliconID", "Chrom", "Start", "End"}
    missing = sorted(required - set(circular_df.columns))
    if missing:
        raise ValueError(f"Missing circular columns: {', '.join(missing)}")

    for amplicon_id, sub in circular_df.groupby("AmpliconID", sort=False):
        sub = sub.reset_index(drop=True).copy()
        n = len(sub)
        for idx, row in sub.iterrows():
            nxt = sub.iloc[(idx + 1) % n]
            c1, p1, clip1 = _segment_exit_breakend(row)
            c2, p2, clip2 = _segment_entry_breakend(nxt)
            rows.append(
                {
                    "AmpliconID": amplicon_id,
                    "JunctionID": junction_id(c1, p1, clip1, c2, p2, clip2),
                    "Chrom1": c1,
                    "Coord1": p1,
                    "Clip1": clip1,
                    "Chrom2": c2,
                    "Coord2": p2,
                    "Clip2": clip2,
                }
            )
    # Keep a stable schema when no circles are inferred.  Downstream matrix
    # construction still needs JunctionID/AmpliconID to write valid empty
    # single-cell outputs instead of failing with KeyError.
    return pd.DataFrame(
        rows,
        columns=[
            "AmpliconID", "JunctionID",
            "Chrom1", "Coord1", "Clip1",
            "Chrom2", "Coord2", "Clip2",
        ],
    )


def same_breakend(a_chrom, a_coord, a_clip, b_chrom, b_coord, b_clip, window=2):
    return str(a_chrom) == str(b_chrom) and str(a_clip) == str(b_clip) and abs(int(a_coord) - int(b_coord)) <= window


def same_junction(bp_row, junction_row, window=2):
    direct = same_breakend(
        bp_row["Chrom1"], bp_row["Coord1"], bp_row["Clip1"],
        junction_row["Chrom1"], junction_row["Coord1"], junction_row["Clip1"],
        window,
    ) and same_breakend(
        bp_row["Chrom2"], bp_row["Coord2"], bp_row["Clip2"],
        junction_row["Chrom2"], junction_row["Coord2"], junction_row["Clip2"],
        window,
    )
    reverse = same_breakend(
        bp_row["Chrom1"], bp_row["Coord1"], bp_row["Clip1"],
        junction_row["Chrom2"], junction_row["Coord2"], junction_row["Clip2"],
        window,
    ) and same_breakend(
        bp_row["Chrom2"], bp_row["Coord2"], bp_row["Clip2"],
        junction_row["Chrom1"], junction_row["Coord1"], junction_row["Clip1"],
        window,
    )
    return direct or reverse


def build_support_long(bp_pair_df, junctions, window=2, cell_barcodes=None):
    barcode_col = "Readbarcode" if "Readbarcode" in bp_pair_df.columns else "Readsbarcode"
    if barcode_col not in bp_pair_df.columns:
        raise ValueError("Breakpoint-pair table must contain Readbarcode or Readsbarcode.")

    rows = []
    for _, bp in bp_pair_df.iterrows():
        barcodes = clean_barcode_list(bp.get(barcode_col))
        if cell_barcodes is not None:
            barcodes = [b for b in barcodes if b in cell_barcodes]
        if not barcodes:
            continue
        for _, junction in junctions.drop_duplicates(["AmpliconID", "JunctionID"]).iterrows():
            if not same_junction(bp, junction, window=window):
                continue
            for barcode, count in Counter(barcodes).items():
                rows.append(
                    {
                        "Barcode": barcode,
                        "AmpliconID": junction["AmpliconID"],
                        "JunctionID": junction["JunctionID"],
                        "support_reads": count,
                    }
                )
    if not rows:
        return pd.DataFrame(columns=["Barcode", "AmpliconID", "JunctionID", "support_reads"])
    out = pd.DataFrame(rows)
    return (
        out.groupby(["Barcode", "AmpliconID", "JunctionID"], sort=False)["support_reads"]
        .sum()
        .reset_index()
    )


def read_spans_coord(read, coord):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.reference_start is None or read.reference_end is None:
        return False
    coord0 = int(coord)
    return read.reference_start <= coord0 < read.reference_end


def read_has_clip_near_coord(read, coord, tolerance=5):
    if not read.cigartuples:
        return False
    coord = int(coord)
    first_op = read.cigartuples[0][0]
    last_op = read.cigartuples[-1][0]
    clipped_ops = {4, 5}
    if first_op in clipped_ops:
        left_clip_coord = int(read.reference_start)
        if abs(left_clip_coord - coord) <= tolerance:
            return True
    if last_op in clipped_ops and read.reference_end is not None:
        right_clip_coord = int(read.reference_end) - 1
        if abs(right_clip_coord - coord) <= tolerance:
            return True
    return False


def is_normal_alignment(read, coord, mapq=20, clip_tolerance=5):
    if read.mapping_quality < mapq:
        return False
    if read.has_tag("SA"):
        return False
    if read_has_clip_near_coord(read, coord, tolerance=clip_tolerance):
        return False
    return read_spans_coord(read, coord)


def count_normal_breakend_barcodes(bam, chrom, coord, search_window=100, mapq=20, barcode_tags=BARCODE_TAGS, clip_tolerance=5, cell_barcodes=None):
    counts = Counter()
    start = max(0, int(coord) - search_window - 1)
    end = int(coord) + search_window
    for read in bam.fetch(str(chrom), start, end):
        if not is_normal_alignment(read, coord, mapq=mapq, clip_tolerance=clip_tolerance):
            continue
        barcode = read_barcode(read, barcode_tags)
        if barcode and (cell_barcodes is None or barcode in cell_barcodes):
            counts[barcode] += 1
    return counts


def build_normal_long(bam_path, junctions, search_window=100, mapq=20, barcode_tags=BARCODE_TAGS, mode="sum", clip_tolerance=5, cell_barcodes=None):
    if pysam is None:
        raise ImportError("pysam is required to build normal-alignment matrices.")
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for _, junction in junctions.drop_duplicates(["AmpliconID", "JunctionID"]).iterrows():
            left = count_normal_breakend_barcodes(
                bam, junction["Chrom1"], junction["Coord1"],
                search_window=search_window, mapq=mapq, barcode_tags=barcode_tags, clip_tolerance=clip_tolerance,
                cell_barcodes=cell_barcodes,
            )
            right = count_normal_breakend_barcodes(
                bam, junction["Chrom2"], junction["Coord2"],
                search_window=search_window, mapq=mapq, barcode_tags=barcode_tags, clip_tolerance=clip_tolerance,
                cell_barcodes=cell_barcodes,
            )
            barcodes = set(left) | set(right)
            for barcode in sorted(barcodes):
                left_count = left.get(barcode, 0)
                right_count = right.get(barcode, 0)
                if mode == "min":
                    normal_reads = min(left_count, right_count)
                else:
                    normal_reads = left_count + right_count
                if normal_reads == 0:
                    continue
                rows.append(
                    {
                        "Barcode": barcode,
                        "AmpliconID": junction["AmpliconID"],
                        "JunctionID": junction["JunctionID"],
                        "normal_reads": normal_reads,
                        "normal_reads_left": left_count,
                        "normal_reads_right": right_count,
                    }
                )
    if not rows:
        return pd.DataFrame(columns=["Barcode", "AmpliconID", "JunctionID", "normal_reads", "normal_reads_left", "normal_reads_right"])
    return pd.DataFrame(rows)


def read_table_or_df(data):
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.read_csv(data, sep="\t")


# ============================================================================ #
# three cell-by-row single-cell matrices
#   (1) cell x ecDNA    binary (ecDNA present / absent)
#   (2) cell x junction  supporting read count
#   (3) cell x junction  summed normal-alignment reads at the two breakends
#       (same rows/cols as (2))
# ============================================================================ #
def _cell_by_key_matrix(long_df, cells, key_ids, value_col):
    """
    Pivot a long table (Barcode / JunctionID / value_col) into a cell x junction
    matrix. Rows fixed to `cells`, columns fixed to `key_ids` (including
    zero-signal junctions), missing filled with 0. Returns a DataFrame whose
    first column is 'Barcode'.
    """
    mat = pd.DataFrame(0, index=pd.Index(cells, name="Barcode"), columns=key_ids)
    if long_df is not None and not long_df.empty:
        # Each shared JunctionID has the same count repeated for each circle.
        # Take it once instead of summing across circle memberships.
        piv = long_df.pivot_table(
            index="Barcode", columns="JunctionID",
            values=value_col, aggfunc="max", fill_value=0,
        )
        mat = piv.reindex(index=cells, columns=key_ids, fill_value=0)
    mat = mat.astype(int)
    mat.index.name = "Barcode"
    return mat.reset_index()


def _cell_by_ecdna_binary(support_long, cells, amplicon_ids, min_support):
    """
    cell x ecDNA binary matrix: 1 if a cell has >= min_support supporting reads on
    ANY junction of that ecDNA, else 0. Rows fixed to `cells`, columns to
    `amplicon_ids`. Returns a DataFrame whose first column is 'Barcode'.
    """
    mat = pd.DataFrame(0, index=pd.Index(cells, name="Barcode"), columns=amplicon_ids)
    if support_long is not None and not support_long.empty:
        amp = (
            support_long.groupby(["Barcode", "AmpliconID"])["support_reads"]
            .max()
            .reset_index()
        )
        amp["pos"] = (amp["support_reads"] >= min_support).astype(int)
        piv = amp.pivot_table(
            index="Barcode", columns="AmpliconID",
            values="pos", aggfunc="max", fill_value=0,
        )
        mat = piv.reindex(index=cells, columns=amplicon_ids, fill_value=0)
    mat = mat.astype(int)
    mat.index.name = "Barcode"
    return mat.reset_index()


def build_cell_matrices(bp_pair, circular, outprefix=None, bam_path=None,
                        window=2, normal_window=100, mapq=20,
                        barcode_tags=BARCODE_TAGS, normal_mode="sum",
                        clip_tolerance=5, min_support=1, cell_barcodes=None,
                        compute_normal=True):
    """
    Build 3 cell-by-row single-cell matrices and (when outprefix is given) write them:
      (1) <prefix>_cell_by_ecDNA_binary.tsv       : cell x ecDNA, 1/0 (cell has ecDNA)
      (2) <prefix>_cell_by_breakpoint_support.tsv : cell x junction, = supporting reads
      (3) <prefix>_cell_by_breakpoint_normal.tsv  : cell x junction, = summed normal
                                                    reads at the two breakends

    - Rows = supplied cell_barcodes, including zero-signal cells; otherwise
      every barcode seen in support or normal.
    - Columns = full set of junctions / ecDNAs inferred from the circle structure
      (including zero-signal ones, so matrix shape stays stable).
    - (2) and (3) share rows/columns (same junctions): one is the circle-junction
      support evidence, the other the reference (normal) evidence.
    - (1) "present" criterion: any junction of that ecDNA has support >= min_support.
    Returns (binary_df, support_df, normal_df). With compute_normal=False,
    skip normal BAM reads and matrix allocation; normal_df is None and no
    normal file is written. Cells then come from support or the whitelist.
    """
    bp_pair = clean_barcode_columns(read_table_or_df(bp_pair))
    circular = read_table_or_df(circular)
    junctions = infer_circular_junctions(circular)

    if cell_barcodes is not None:
        cell_barcodes = set(cell_barcodes)
    support_long = build_support_long(bp_pair, junctions, window=window,
                                      cell_barcodes=cell_barcodes)
    if bam_path and compute_normal:
        normal_long = build_normal_long(
            bam_path, junctions, search_window=normal_window, mapq=mapq,
            barcode_tags=barcode_tags, mode=normal_mode, clip_tolerance=clip_tolerance,
            cell_barcodes=cell_barcodes,
        )
    else:
        normal_long = pd.DataFrame(
            columns=["Barcode", "AmpliconID", "JunctionID", "normal_reads"]
        )

    # rows: cells with any signal (barcodes appearing in support or normal)
    cells = sorted(
        set(support_long["Barcode"].tolist()) | set(normal_long["Barcode"].tolist())
    )
    if cell_barcodes is not None:
        cells = sorted(set(cell_barcodes))
    # columns: full junction / ecDNA set (keep order of appearance in the circle)
    junction_ids = list(dict.fromkeys(junctions["JunctionID"].tolist()))
    amplicon_ids = list(dict.fromkeys(junctions["AmpliconID"].tolist()))

    support_df = _cell_by_key_matrix(support_long, cells, junction_ids, "support_reads")
    normal_df = (_cell_by_key_matrix(normal_long, cells, junction_ids, "normal_reads")
                 if compute_normal else None)
    binary_df = _cell_by_ecdna_binary(support_long, cells, amplicon_ids, min_support)

    if outprefix is not None:
        outprefix = Path(outprefix)
        outprefix.parent.mkdir(parents=True, exist_ok=True)
        binary_df.to_csv(f"{outprefix}_cell_by_ecDNA_binary.tsv", sep="\t", index=False)
        support_df.to_csv(f"{outprefix}_cell_by_breakpoint_support.tsv", sep="\t", index=False)
        if normal_df is not None:
            normal_df.to_csv(f"{outprefix}_cell_by_breakpoint_normal.tsv", sep="\t", index=False)

    return binary_df, support_df, normal_df




def write_cell_matrices_if_single_cell(bp_pair, circular, outprefix, bam_path=None,
                                       window=2, normal_window=100, mapq=20,
                                       barcode_tags=BARCODE_TAGS, normal_mode="sum",
                                       clip_tolerance=5, min_support=1):
    """Build the 3 cell matrices only when the bp-pair table has single-cell
    barcodes, else return None."""
    bp_pair_df = clean_barcode_columns(read_table_or_df(bp_pair))
    if not has_single_cell_barcodes(bp_pair_df):
        return None
    return build_cell_matrices(
        bp_pair_df, circular, outprefix, bam_path=bam_path, window=window,
        normal_window=normal_window, mapq=mapq, barcode_tags=barcode_tags,
        normal_mode=normal_mode, clip_tolerance=clip_tolerance, min_support=min_support,
    )


def write_cell_assignments_if_single_cell(bp_pair, circular, outprefix,
                                         bam_path=None, fragments_path=None,
                                         cell_barcodes=None, derived_circular=None,
                                         raw_bp_pair=None, all_min_count=5, **kwargs):
    """Write and return (cell/ecDNA binary, ALL raw-pair support, ALL normal).

    ALL matrices retain all whitelist cells (or raw/normal observed cells),
    alongside binary assignment. raw_bp_pair is required.
    ALL uses Count>=all_min_count and the instability reference normal filters.
    Write breakpoint_annotation_ALL/ecDNA tables keyed by matrix JunctionID.

    Binary assignment requires fragments and retains positive and negative
    cells; without a whitelist it uses all fragment barcodes.
    Derived circles contribute junction columns only.
    Bulk inputs without barcodes or explicit single-cell inputs are skipped.
    Also write <prefix>_cell_assignment_evidence_Q1_no_pair.tsv, retaining
    per-ring Q1 thresholds and reasons for the final binary decisions.
    No fixed-cutoff baseline is computed. The filename is retained for
    compatibility; effective_rule_version identifies the replacement rule.
    """
    bp_pair = clean_barcode_columns(read_table_or_df(bp_pair))
    if raw_bp_pair is None:
        raise ValueError("ALL matrices require raw_bp_pair (the pre-graph raw pair table)")
    raw_bp_pair = clean_barcode_columns(read_table_or_df(raw_bp_pair))
    if (not has_single_cell_barcodes(raw_bp_pair) and fragments_path is None
            and cell_barcodes is None):
        return None
    try:
        from .sc_assignment import assign_cells
    except ImportError:
        from sc_assignment import assign_cells
    # Q1 uses support, weak evidence and fragments; ALL normal is built below.
    kwargs['compute_normal'] = False
    binary, support, normal, evidence = assign_cells(
        bp_pair, circular, bam_path=bam_path, fragments_path=fragments_path,
        cell_barcodes=cell_barcodes, derived_circular=derived_circular,
        return_evidence=True, **kwargs)
    try:
        from .sc_breakpoint import build_all_breakpoint_matrices
    except ImportError:
        from sc_breakpoint import build_all_breakpoint_matrices
    circles = read_table_or_df(circular)
    derived = derived_circular
    if derived is None:
        metadata = circles.attrs.get("circle_dedup", {})
        derived = pd.DataFrame(metadata.get("derived_segments", []),
                               columns=metadata.get("segment_columns", circles.columns))
    else:
        derived = read_table_or_df(derived)
    all_circles = pd.concat([circles, derived], ignore_index=True)
    support, normal, _, _ = build_all_breakpoint_matrices(
        raw_bp_pair, all_circles, bam_path, outprefix=outprefix,
        cell_barcodes=cell_barcodes, min_count=all_min_count, mapq=20)
    matrices = (binary, support, normal)
    prefix = Path(outprefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for kind, table in zip(("ecDNA_binary", "breakpoint_support", "breakpoint_normal"), matrices):
        table.to_csv(f"{prefix}_cell_by_{kind}.tsv", sep="\t", index=False)
    evidence.to_csv(f"{prefix}_cell_assignment_evidence_Q1_no_pair.tsv", sep="\t", index=False)
    summary = binary.attrs["assignment_summary"]
    Path(f"{prefix}_cell_assignment_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Cell assignment: {summary['assignment_mode']}; candidates={summary['candidate_cells']}; "
          f"returned={summary['returned_cells']}; positive_only={summary['positive_only']}")
    return matrices
