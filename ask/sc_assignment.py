"""Single-cell assignment for retained ASK circles (no additional Score filter).

Weak/region evidence helpers adapted from ask_sc_evidence_ask.py (20260825).
Input circular rows must be the retained representatives in structural order.
This module performs no output-file writes and never uses simulation truth.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import gzip
import math
import numpy as np
import pandas as pd
import pysam

BARCODE_TAGS = ("CB", "CR", "BX", "BC")
def read_barcode(read: pysam.AlignedSegment) -> str | None:
    for tag in BARCODE_TAGS:
        if read.has_tag(tag):
            value = str(read.get_tag(tag)).strip()
            return value or None
    return None

def correct_clip(read: pysam.AlignedSegment, coordinate: int, clip: str, tolerance: int) -> bool:
    if not read.cigartuples:
        return False
    clipped_operations = {4, 5}
    if clip == "L":
        return (read.cigartuples[0][0] in clipped_operations and abs(int(read.reference_start) - coordinate) <= tolerance)
    return (read.cigartuples[-1][0] in clipped_operations
            and read.reference_end is not None
            and abs(int(read.reference_end) - coordinate) <= tolerance)

def sa_hits_partner(read: pysam.AlignedSegment, chromosome: str, coordinate: int, window: int) -> bool:
    if not read.has_tag("SA"):
        return False
    for item in str(read.get_tag("SA")).split(";"):
        fields = item.split(",")
        if len(fields) >= 2 and fields[0] == chromosome:
            if abs((int(fields[1]) - 1) - coordinate) <= window:
                return True
    return False

def mate_hits_partner(read: pysam.AlignedSegment, bam: pysam.AlignmentFile,
                      chromosome: str, coordinate: int, window: int) -> bool:
    return (
        not read.mate_is_unmapped
        and read.next_reference_id >= 0
        and bam.get_reference_name(read.next_reference_id) == chromosome
        and abs(int(read.next_reference_start) - coordinate) <= window
    )

def scan_weak_evidence(bam_path: Path, junctions: pd.DataFrame,
                       fetch_window: int, partner_window: int,
                       clip_tolerance: int, min_mapq: int,
                       selected_cells: set[str]) -> pd.DataFrame:
    """Scan BAM only around the junctions of retained high-quality circles."""
    evidence = {}
    weights = {"split_sa": 1.0, "discordant_pair": 0.7, "single_end_clip": 0.3}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for junction in junctions.itertuples(index=False):
            sides = [
                (junction.Chrom1, int(junction.Coord1), junction.Clip1,
                 junction.Chrom2, int(junction.Coord2)),
                (junction.Chrom2, int(junction.Coord2), junction.Clip2,
                 junction.Chrom1, int(junction.Coord1)),
            ]
            for chrom, coord, clip, partner_chrom, partner_coord in sides:
                try:
                    reads = bam.fetch(
                        chrom, max(0, coord - fetch_window), coord + fetch_window + 1)
                except (ValueError, KeyError):
                    continue
                for read in reads:
                    if (read.is_unmapped or read.is_secondary or read.mapping_quality < min_mapq):
                        continue
                    barcode = read_barcode(read)
                    if barcode not in selected_cells or (read.has_tag("NH") and int(read.get_tag("NH")) > 1):
                        continue
                    clipped = correct_clip(read, coord, clip, clip_tolerance)
                    if clipped and sa_hits_partner(read, partner_chrom, partner_coord, partner_window):
                        evidence_type = "split_sa"
                    elif mate_hits_partner(read, bam, partner_chrom, partner_coord, partner_window):
                        evidence_type = "discordant_pair"
                    elif clipped:
                        evidence_type = "single_end_clip"
                    else:
                        continue
                    key = (barcode, str(junction.AmpliconID), str(junction.JunctionID), str(read.query_name))
                    previous = evidence.get(key)
                    if previous is None or weights[evidence_type] > weights[previous]:
                        evidence[key] = evidence_type
    rows = [
        {"Barcode": barcode, "AmpliconID": amp, "JunctionID": junction,
         "ReadName": read_name, "evidence_type": evidence_type,
         "weight": weights[evidence_type]}
        for (barcode, amp, junction, read_name), evidence_type in evidence.items()
    ]
    return pd.DataFrame(rows, columns=[ "Barcode", "AmpliconID", "JunctionID", "ReadName",
        "evidence_type", "weight" ])

def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()

def make_flanks(circular: pd.DataFrame) -> dict[str, list[tuple[str, int, int]]]:
    ring_intervals = defaultdict(list)
    for row in circular.itertuples(index=False):
        ring_intervals[str(row.Chrom)].append((int(row.Start), int(row.End)))
    output = defaultdict(list)
    for row in circular.itertuples(index=False):
        chrom, start, end = str(row.Chrom), int(row.Start), int(row.End)
        length = max(1, end - start)
        for left, right in ((max(0, start - length), start), (end, end + length)):
            overlaps_circle = any(
                max(left, other_left) < min(right, other_right)
                for other_left, other_right in ring_intervals[chrom])
            if not overlaps_circle and right > left:
                output[str(row.AmpliconID)].append((chrom, left, right))
    return output

def fragment_barcodes(fragment_path: Path) -> set[str]:
    """Collect every barcode in the fourth column of a fragments file."""
    barcodes = set()
    with open_text(fragment_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3].strip():
                barcodes.add(fields[3].strip())
    if not barcodes:
        raise ValueError(f"No barcodes found in fragments: {fragment_path}")
    return barcodes


def fragment_features(fragment_path: Path, circular: pd.DataFrame,
                      selected_cells: set[str]) -> tuple[Counter, dict, dict, dict, dict]:
    segments = defaultdict(list)
    segment_keys = defaultdict(list)
    for index, row in circular.reset_index(drop=True).iterrows():
        key = f"{row.AmpliconID}:seg{index}"
        segments[str(row.Chrom)].append(
            (int(row.Start), int(row.End), str(row.AmpliconID), key))
        segment_keys[str(row.AmpliconID)].append(key)
    flanks_by_amp = make_flanks(circular)
    flanks = defaultdict(list)
    for amp, intervals in flanks_by_amp.items():
        for chrom, start, end in intervals:
            flanks[chrom].append((start, end, amp))

    total = Counter()
    segment_counts = defaultdict(Counter)
    flank_counts = defaultdict(Counter)
    with open_text(fragment_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom, start_text, end_text, barcode = fields[:4]
            if barcode not in selected_cells:
                continue
            try:
                start, end = int(start_text), int(end_text)
                count = int(fields[4]) if len(fields) > 4 else 1
            except ValueError:
                continue
            total[barcode] += count
            for left, right, amp, key in segments.get(chrom, []):
                if start < right and end > left:
                    segment_counts[(barcode, amp)][key] += count
            for left, right, amp in flanks.get(chrom, []):
                if start < right and end > left:
                    flank_counts[(barcode, amp)][f"{chrom}:{left}-{right}"] += count
    return total, segment_counts, flank_counts, dict(segment_keys), flanks_by_amp

def make_features(cells: list[str], amps: list[str], direct: Counter,
                  distinct: dict, weak: pd.DataFrame, total: Counter,
                  segment_counts: dict, flank_counts: dict,
                  segment_keys: dict, ring_lengths: dict,
                  flank_lengths: dict, sample_name: str) -> pd.DataFrame:
    weak_score = Counter()
    weak_count = Counter()
    weak_types = defaultdict(Counter)
    for row in weak.itertuples(index=False):
        key = (str(row.Barcode), str(row.AmpliconID))
        weak_score[key] += float(row.weight)
        weak_count[key] += 1
        weak_types[key][str(row.evidence_type)] += 1

    rows = []
    for barcode in cells:
        for amp in amps:
            key = (barcode, amp)
            ordered_keys = segment_keys[amp]
            values = [int(segment_counts.get(key, {}).get(segment_key, 0))
                      for segment_key in ordered_keys]
            flank_values = list(flank_counts.get(key, {}).values())
            segment_count = int(sum(values))
            flank_count = int(sum(flank_values))
            total_fragments = int(total.get(barcode, 0))
            ring_kb = max(float(ring_lengths[amp]), 1.0) / 1000.0
            flank_kb = max(float(flank_lengths[amp]), 1.0) / 1000.0
            segment_rate = (
                1e6 * segment_count / total_fragments / ring_kb
                if total_fragments else 0.0)
            flank_rate = (
                1e6 * flank_count / total_fragments / flank_kb
                if total_fragments else 0.0)
            covered = sum(value > 0 for value in values)
            rows.append({
                "Sample": sample_name, "Barcode": barcode, "AmpliconID": amp,
                "n_junction_reads": int(direct.get(key, 0)),
                "n_distinct_junctions": len(distinct.get(key, set())),
                "n_weak_reads": int(weak_count.get(key, 0)),
                "weak_score": float(weak_score.get(key, 0.0)),
                "n_weak_split_sa": weak_types[key]["split_sa"],
                "n_weak_discordant_pairs": weak_types[key]["discordant_pair"],
                "n_weak_single_clips": weak_types[key]["single_end_clip"],
                "segment_count": segment_count, "flank_count": flank_count,
                "total_fragments": total_fragments,
                "segment_rate_per_million_per_kb": segment_rate,
                "flank_rate_per_million_per_kb": flank_rate,
                "segment_log_enrichment": math.log1p(segment_rate) - math.log1p(flank_rate),
                "n_segments_total": len(values),
                "n_segments_covered": covered,
                "fraction_segments_covered": covered / len(values) if values else 0.0,
            })
    features = pd.DataFrame(rows)

    anchor_profiles = {}
    for amp in amps:
        normalized = []
        for barcode in cells:
            if direct.get((barcode, amp), 0) < 1:
                continue
            vector = np.array([
                segment_counts.get((barcode, amp), {}).get(key, 0)
                for key in segment_keys[amp]
            ], dtype=float)
            if vector.sum() > 0:
                normalized.append(vector / vector.sum())
        if normalized:
            anchor_profiles[amp] = np.mean(normalized, axis=0)

    similarities = []
    for row in features.itertuples(index=False):
        vector = np.array([
            segment_counts.get((str(row.Barcode), str(row.AmpliconID)), {}).get(key, 0)
            for key in segment_keys[str(row.AmpliconID)]
        ], dtype=float)
        anchor = anchor_profiles.get(str(row.AmpliconID))
        if anchor is None or np.linalg.norm(vector) == 0 or np.linalg.norm(anchor) == 0:
            similarities.append(0.0)
        else:
            similarities.append(float(
                np.dot(vector, anchor) / (np.linalg.norm(vector) * np.linalg.norm(anchor))))
    features["coverage_profile_cosine"] = similarities
    return features

def assign_rules(features: pd.DataFrame, weak_score_threshold: float,
                 ring_q1: pd.Series | float, log_enrichment_threshold: float,
                 profile_cosine_threshold: float,
                 segment_fraction_threshold: float) -> pd.DataFrame:
    output = features.copy()
    direct = output["n_junction_reads"].ge(1)
    weak = ~direct & output["weak_score"].ge(weak_score_threshold)
    multi = (
        ~direct & ~weak
        & output["segment_rate_per_million_per_kb"].ge(ring_q1)
        & output["segment_log_enrichment"].ge(log_enrichment_threshold)
        & output["coverage_profile_cosine"].ge(profile_cosine_threshold)
        & output["fraction_segments_covered"].ge(segment_fraction_threshold)
    )
    output["assignment"] = np.select(
        [direct, weak | multi],
        ["high_confidence_positive", "probable_positive"],
        default="negative",
    )
    output["assignment_reason"] = np.select(
        [direct, weak, multi],
        ["final_junction_read", "targeted_weak_breakpoint_evidence",
         "multi_segment_rate_enrichment"],
        default="insufficient_evidence",
    )
    output["ecDNA_binary"] = (direct | weak | multi).astype(int)
    return output


def rescue_ring_q1(evidence, log_enrichment_threshold=2,
                   weak_score_threshold=0.3, profile_cosine_threshold=0.96,
                   segment_fraction_threshold=1):
    """Reclassify with Q1 replacing the fixed coverage cutoff (weak05 rule).

    Despite the retained public function name, this is not additive rescue.
    Keep direct/weak calls and require all coverage gates for the Q1 path.
    Previous coverage-only positives may be removed. Q1 is computed from
    direct-positive cells within each sample and ring; no reference disables
    the coverage path without a fallback to the fixed cutoff.
    """
    output = evidence.copy()
    group_columns = (["Sample"] if "Sample" in output.columns else []) + ["AmpliconID"]
    direct = output["n_junction_reads"].ge(1)
    thresholds = output.loc[direct].groupby(group_columns, dropna=False)[
        "segment_rate_per_million_per_kb"].quantile(0.25)
    if len(group_columns) == 1:
        output["ring_direct_Q1"] = output["AmpliconID"].map(thresholds)
    else:
        keys = pd.MultiIndex.from_frame(output[group_columns])
        output["ring_direct_Q1"] = thresholds.reindex(keys).to_numpy()
    obsolete = [c for c in output if c.startswith("original_") or c in
                ("ring_Q1_rescue_added", "ring_Q1_replacement_removed")]
    output = output.drop(columns=obsolete)
    updated = assign_rules(
        output, weak_score_threshold, output["ring_direct_Q1"],
        log_enrichment_threshold, profile_cosine_threshold,
        segment_fraction_threshold)
    coverage = updated["assignment_reason"].eq("multi_segment_rate_enrichment")
    updated.loc[coverage, "assignment_reason"] = "per_ring_Q1_coverage_replacement"
    updated["ring_Q1_log_enrichment_threshold"] = log_enrichment_threshold
    updated["effective_rule_version"] = "per_ring_Q1_only_weak05_20261110"
    return updated


def assign_cells(bp_pair, circular, bam_path=None, fragments_path=None,
                 window=2, normal_window=100, mapq=20,
                 cell_barcodes=None, derived_circular=None,
                 weak_score_threshold=0.3,
                 log_enrichment_threshold=2, profile_cosine_threshold=0.96,
                 segment_fraction_threshold=1,
                 return_evidence=False, compute_normal=True):
    """Return binary, support and normal matrices for single-cell assignment.

    Fragments are required. Use all fragment barcodes unless a whitelist is
    supplied, and return positive and negative cells with Q1 assignment.
    Binary columns are retained representatives;
    junction columns include representatives and derived circles, counted once.
    With return_evidence=True, append the per-cell/ring decision table to the
    three matrices. With fragments, Q1 replaces the fixed coverage cutoff
    and retains enrichment, cosine, and covered-fraction gates (weak05 rule).
    No fixed coverage threshold or baseline classification is computed.
    compute_normal=False skips the legacy normal scan and returns None for
    normal; it does not change Q1 decisions or the candidate-cell universe.
    """
    try:
        from . import sc_matrix as sc
    except ImportError:
        import sc_matrix as sc
    if fragments_path is not None and not Path(fragments_path).is_file():
        raise FileNotFoundError(fragments_path)
    if fragments_path is None:
        try:
            from .sc_fragments import report_missing_fragments
        except ImportError:
            from sc_fragments import report_missing_fragments
        report_missing_fragments(bam_path)
        raise ValueError("Fragments are required for single-cell assignment; provide fragments_path.")
    bp_pair = sc.clean_barcode_columns(sc.read_table_or_df(bp_pair))
    if cell_barcodes is not None:
        allowed = {str(b).strip() for b in cell_barcodes if str(b).strip()}
        if not allowed:
            raise ValueError("--sc-barcodes must contain at least one barcode")
        candidate_source = "whitelist"
    else:
        allowed = fragment_barcodes(Path(fragments_path))
        candidate_source = "fragments"
    mode = "Q1"
    metadata = dict(assignment_mode=mode, positive_only=False,
                    candidate_source=candidate_source,
                    candidate_cells=len(allowed), coverage_available=True)
    if bam_path is None:
        raise ValueError("Cell assignment requires an indexed BAM")
    circular = sc.read_table_or_df(circular).copy()
    junctions = sc.infer_circular_junctions(circular)
    amps = list(dict.fromkeys(circular["AmpliconID"].astype(str)))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if not bam.has_index():
            raise ValueError("Cell assignment requires an indexed BAM")
    if derived_circular is None:
        circle_metadata = circular.attrs.get("circle_dedup", {})
        derived_circular = pd.DataFrame(circle_metadata.get("derived_segments", []),
                                        columns=circle_metadata.get("segment_columns", circular.columns))
    else:
        derived_circular = sc.read_table_or_df(derived_circular)
    all_circles = pd.concat([circular, derived_circular], ignore_index=True)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        missing = set(all_circles.Chrom.astype(str)) - set(bam.references)
        if missing:
            raise ValueError(f"Circle chromosomes absent from BAM: {sorted(missing)}")
    _, support, normal = sc.build_cell_matrices(
        bp_pair, all_circles, bam_path=bam_path, window=window,
        normal_window=normal_window, mapq=mapq, cell_barcodes=allowed,
        compute_normal=compute_normal)
    support = support.set_index("Barcode")
    if normal is not None:
        normal = normal.set_index("Barcode")
    cells = set(allowed)
    if not amps or not cells:
        output_cells = sorted(cells)
        binary = pd.DataFrame(0, index=pd.Index(output_cells, name="Barcode"), columns=amps).reset_index()
        binary.attrs["assignment_summary"] = dict(metadata, returned_cells=len(output_cells), positive_cells=0)
        matrices = (binary, support.reindex(output_cells).reset_index(),
                    normal.reindex(output_cells).reset_index() if normal is not None else None)
        if return_evidence:
            evidence = pd.DataFrame(columns=[
                "Barcode", "AmpliconID", "ecDNA_binary", "assignment",
                "assignment_reason", "ring_direct_Q1",
                "effective_rule_version",
            ])
            return (*matrices, evidence)
        return matrices
    direct, distinct = Counter(), defaultdict(set)
    for amp in amps:
        ids = list(dict.fromkeys(junctions.loc[
            junctions.AmpliconID.astype(str) == amp, "JunctionID"]))
        for barcode, row in support[ids].iterrows():
            if row.sum() > 0:
                direct[(barcode, amp)] = int(row.sum())
                distinct[(barcode, amp)] = set(row.index[row > 0])
    weak = scan_weak_evidence(Path(bam_path), junctions, 100, 250, 10, 10, allowed)
    total, counts, flank_counts, keys, flanks = fragment_features(Path(fragments_path), circular, allowed)
    lengths = {str(a): int((b.End - b.Start).sum())
               for a, b in circular.groupby("AmpliconID")}
    flank_lengths = {a: sum(e-s for _, s, e in flanks.get(a, [])) or 1 for a in amps}
    features = make_features(sorted(cells), amps, direct, distinct, weak,
                             total, counts, flank_counts, keys, lengths,
                             flank_lengths, "")
    evidence = rescue_ring_q1(
        features, log_enrichment_threshold,
        weak_score_threshold=weak_score_threshold,
        profile_cosine_threshold=profile_cosine_threshold,
        segment_fraction_threshold=segment_fraction_threshold)
    evidence["effective_rule_version"] = "per_ring_Q1_only_weak05_20261110"
    for key, value in metadata.items():
        evidence[key] = value
    decision_evidence = evidence.copy()
    evidence = evidence.rename(columns={"AmpliconID": "ecDNA_ID"})
    cells = sorted(cells)
    support = support.reindex(cells, fill_value=0)
    if normal is not None:
        normal = normal.reindex(cells, fill_value=0)
    binary = evidence.pivot(index="Barcode", columns="ecDNA_ID",
                            values="ecDNA_binary").reindex(index=cells, columns=amps, fill_value=0)
    binary = binary.fillna(0).astype(int)
    binary.columns.name = None
    positive_cells = int(binary.any(axis=1).sum())
    for frame in (binary, support, normal):
        if frame is not None:
            frame.index.name = "Barcode"
    matrices = tuple(t.reset_index() if t is not None else None
                     for t in (binary, support, normal))
    matrices[0].attrs["assignment_summary"] = dict(metadata, returned_cells=len(binary), positive_cells=positive_cells)
    return (*matrices, decision_evidence) if return_evidence else matrices
