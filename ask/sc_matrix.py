import ast
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


def infer_circular_junctions(circular_df):
    rows = []
    required = {"AmpliconID", "Chrom", "Start", "End"}
    missing = sorted(required - set(circular_df.columns))
    if missing:
        raise ValueError(f"Missing circular columns: {', '.join(missing)}")

    for amplicon_id, sub in circular_df.groupby("AmpliconID", sort=False):
        sub = sub.reset_index(drop=True).copy()
        for idx, row in sub.iterrows():
            if len(sub) == 1:
                rows.append(
                    {
                        "AmpliconID": amplicon_id,
                        "JunctionID": junction_id(row["Chrom"], row["Start"], "L", row["Chrom"], row["End"], "R"),
                        "Chrom1": row["Chrom"],
                        "Coord1": int(row["Start"]),
                        "Clip1": "L",
                        "Chrom2": row["Chrom"],
                        "Coord2": int(row["End"]),
                        "Clip2": "R",
                    }
                )
                continue
            nxt = sub.iloc[(idx + 1) % len(sub)]
            rows.append(
                {
                    "AmpliconID": amplicon_id,
                    "JunctionID": junction_id(row["Chrom"], row["End"], "R", nxt["Chrom"], nxt["Start"], "L"),
                    "Chrom1": row["Chrom"],
                    "Coord1": int(row["End"]),
                    "Clip1": "R",
                    "Chrom2": nxt["Chrom"],
                    "Coord2": int(nxt["Start"]),
                    "Clip2": "L",
                }
            )
    return pd.DataFrame(rows)


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


def build_support_long(bp_pair_df, junctions, window=2):
    barcode_col = "Readbarcode" if "Readbarcode" in bp_pair_df.columns else "Readsbarcode"
    if barcode_col not in bp_pair_df.columns:
        raise ValueError("Breakpoint-pair table must contain Readbarcode or Readsbarcode.")

    rows = []
    for _, bp in bp_pair_df.iterrows():
        barcodes = clean_barcode_list(bp.get(barcode_col))
        if not barcodes:
            continue
        for _, junction in junctions.iterrows():
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


def long_to_matrix(long_df, value_col):
    if long_df.empty:
        return pd.DataFrame()
    matrix = (
        long_df.pivot_table(index="JunctionID", columns="Barcode", values=value_col, aggfunc="sum", fill_value=0)
        .reset_index()
        .rename_axis(None, axis=1)
    )
    amplicons = long_df[["JunctionID", "AmpliconID"]].drop_duplicates("JunctionID")
    return amplicons.merge(matrix, on="JunctionID", how="right")


def align_matrices(junctions, support_matrix, normal_matrix):
    key_cols = ["JunctionID", "AmpliconID"]
    template = junctions[key_cols].drop_duplicates("JunctionID").reset_index(drop=True)
    barcode_cols = sorted((set(support_matrix.columns) | set(normal_matrix.columns)) - set(key_cols))

    def align_one(matrix):
        out = template.copy()
        if not matrix.empty:
            out = out.merge(matrix, on=key_cols, how="left")
        for barcode in barcode_cols:
            if barcode not in out.columns:
                out[barcode] = 0
        if barcode_cols:
            out[barcode_cols] = out[barcode_cols].fillna(0).astype(int)
        return out[key_cols + barcode_cols]

    return align_one(support_matrix), align_one(normal_matrix)


def read_spans_coord(read, coord_1based):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.reference_start is None or read.reference_end is None:
        return False
    coord0 = int(coord_1based) - 1
    return read.reference_start <= coord0 < read.reference_end


def read_has_clip_near_coord(read, coord_1based, tolerance=5):
    if not read.cigartuples:
        return False
    coord = int(coord_1based)
    first_op = read.cigartuples[0][0]
    last_op = read.cigartuples[-1][0]
    clipped_ops = {4, 5}
    if first_op in clipped_ops:
        left_clip_coord = int(read.reference_start) + 1
        if abs(left_clip_coord - coord) <= tolerance:
            return True
    if last_op in clipped_ops and read.reference_end is not None:
        right_clip_coord = int(read.reference_end)
        if abs(right_clip_coord - coord) <= tolerance:
            return True
    return False


def is_normal_alignment(read, coord_1based, mapq=20, clip_tolerance=5):
    if read.mapping_quality < mapq:
        return False
    if read.has_tag("SA"):
        return False
    if read_has_clip_near_coord(read, coord_1based, tolerance=clip_tolerance):
        return False
    return read_spans_coord(read, coord_1based)


def count_normal_breakend_barcodes(bam, chrom, coord, search_window=100, mapq=20, barcode_tags=BARCODE_TAGS, clip_tolerance=5):
    counts = Counter()
    start = max(0, int(coord) - search_window - 1)
    end = int(coord) + search_window
    for read in bam.fetch(str(chrom), start, end):
        if not is_normal_alignment(read, coord, mapq=mapq, clip_tolerance=clip_tolerance):
            continue
        barcode = read_barcode(read, barcode_tags)
        if barcode:
            counts[barcode] += 1
    return counts


def build_normal_long(bam_path, junctions, search_window=100, mapq=20, barcode_tags=BARCODE_TAGS, mode="sum", clip_tolerance=5):
    if pysam is None:
        raise ImportError("pysam is required to build normal-alignment matrices.")
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for _, junction in junctions.iterrows():
            left = count_normal_breakend_barcodes(
                bam, junction["Chrom1"], junction["Coord1"],
                search_window=search_window, mapq=mapq, barcode_tags=barcode_tags, clip_tolerance=clip_tolerance,
            )
            right = count_normal_breakend_barcodes(
                bam, junction["Chrom2"], junction["Coord2"],
                search_window=search_window, mapq=mapq, barcode_tags=barcode_tags, clip_tolerance=clip_tolerance,
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


def build_sc_matrices(bp_pair, circular, outprefix=None, bam_path=None, window=2, normal_window=100, mapq=20, barcode_tags=BARCODE_TAGS, normal_mode="sum", clip_tolerance=5):
    bp_pair = clean_barcode_columns(read_table_or_df(bp_pair))
    circular = read_table_or_df(circular)

    junctions = infer_circular_junctions(circular)
    support_long = build_support_long(bp_pair, junctions, window=window)
    support_matrix = long_to_matrix(support_long, "support_reads")

    if bam_path:
        normal_long = build_normal_long(
            bam_path,
            junctions,
            search_window=normal_window,
            mapq=mapq,
            barcode_tags=barcode_tags,
            mode=normal_mode,
            clip_tolerance=clip_tolerance,
        )
        normal_matrix = long_to_matrix(normal_long, "normal_reads")
    else:
        normal_matrix = pd.DataFrame()

    support_matrix, normal_matrix = align_matrices(junctions, support_matrix, normal_matrix)
    if outprefix is not None:
        outprefix = Path(outprefix)
        outprefix.parent.mkdir(parents=True, exist_ok=True)
        support_matrix.to_csv(f"{outprefix}_support_matrix.tsv", sep="\t", index=False)
        normal_matrix.to_csv(f"{outprefix}_normal_alignment_matrix.tsv", sep="\t", index=False)

    return support_matrix, normal_matrix


def write_sc_matrices_if_single_cell(bp_pair, circular, outprefix, bam_path=None, window=2, normal_window=100, mapq=20, barcode_tags=BARCODE_TAGS, normal_mode="sum", clip_tolerance=5):
    bp_pair_df = clean_barcode_columns(read_table_or_df(bp_pair))
    if not has_single_cell_barcodes(bp_pair_df):
        return None
    return build_sc_matrices(
        bp_pair_df,
        circular,
        outprefix,
        bam_path=bam_path,
        window=window,
        normal_window=normal_window,
        mapq=mapq,
        barcode_tags=barcode_tags,
        normal_mode=normal_mode,
        clip_tolerance=clip_tolerance,
    )
