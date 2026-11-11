"""Read-level evidence normalization and junction-level support merging."""

import numpy as np
import pandas as pd

import covbam


BP_PAIR_COLUMNS = [
    "Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", "Count",
    "offset", "Seq", "ReadID", "Readbarcode", "HasPEEvidence",
]


def _evidence_flag(value):
    """Return a stable bool for flags loaded either in memory or from TSV."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", "y"}
    return bool(value)


def _evidence_list(value):
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return [value.item()]
        return value.tolist()
    if isinstance(value, (list, tuple, set, pd.Series)):
        return list(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [value]


def _read_id_array(value):
    """Normalize uint64 arrays and upgrade legacy QNAME lists when needed."""
    values = _evidence_list(value)
    if not values:
        return np.empty(0, dtype=np.uint64)

    if all(isinstance(item, (int, np.integer)) for item in values):
        return np.asarray(values, dtype=np.uint64)

    read_ids, _ = covbam.encode_read_evidence(values)
    return read_ids


def _upgrade_read_evidence(dataframe):
    """Ensure one compact ReadID array and aligned barcode list per evidence row."""
    dataframe = dataframe.copy()
    if 'ReadID' not in dataframe.columns and 'Readname' in dataframe.columns:
        if 'Readbarcode' not in dataframe.columns:
            dataframe['Readbarcode'] = [[] for _ in range(dataframe.shape[0])]
        encoded = [
            covbam.encode_read_evidence(readnames, readbarcodes)
            for readnames, readbarcodes in zip(dataframe['Readname'], dataframe['Readbarcode'])
        ]
        dataframe['ReadID'] = [item[0] for item in encoded]
        dataframe['Readbarcode'] = [item[1] for item in encoded]
    elif 'ReadID' not in dataframe.columns:
        dataframe['ReadID'] = [np.empty(0, dtype=np.uint64) for _ in range(dataframe.shape[0])]
    else:
        dataframe['ReadID'] = dataframe['ReadID'].map(_read_id_array)

    if 'Readbarcode' not in dataframe.columns:
        dataframe['Readbarcode'] = [[] for _ in range(dataframe.shape[0])]

    return dataframe.drop(columns=['Readname'], errors='ignore')


def _canonical_junction_key(row):
    end1 = (str(row['Chrom1']), int(row['Coord1']), str(row['Clip1']))
    end2 = (str(row['Chrom2']), int(row['Coord2']), str(row['Clip2']))
    return (end1, end2) if end1 <= end2 else (end2, end1)


def edges_to_bp_pair_df(edges):
    """Convert discordant edges to breakpoint-pair rows with compact ReadID arrays."""
    rows = []
    for edge in edges:
        if len(edge) < 3:
            continue

        (chromA, posA, strA), (chromB, posB, strB), count = edge[:3]
        read_ids = _read_id_array(edge[3] if len(edge) > 3 else [])
        readbarcodes = _evidence_list(edge[4] if len(edge) > 4 else [])
        readbarcodes = [str(barcode) if barcode is not None and str(barcode) not in {'', 'None', 'nan'} else None for barcode in readbarcodes]
        # BARCODE: [] means this WGS edge has no barcode; do not expand it to N None values.
        if readbarcodes and len(readbarcodes) < read_ids.size:
            readbarcodes.extend([None] * (read_ids.size - len(readbarcodes)))
        elif readbarcodes and len(readbarcodes) > read_ids.size:
            readbarcodes = readbarcodes[:read_ids.size]

        clipA = 'R' if strA == '+' else 'L'
        clipB = 'R' if strB == '+' else 'L'
        support = int(count)
        rows.append([
            chromA, posA, clipA, chromB, posB, clipB, support, 0, 'DS',
            read_ids, readbarcodes, True,
        ])

    return pd.DataFrame(rows, columns=BP_PAIR_COLUMNS)


def merge_bp_pair_evidence(*dataframes):
    """Merge evidence, using earlier dataframes as higher-priority sources.

    Read evidence is combined across every matching row and deduplicated by the
    stable QNAME-derived ``ReadID``.  Therefore ``Count`` is the number of
    unique supporting read pairs/molecules for rows with tracked read evidence,
    rather than the sum of split/PE/DS alignment occurrences.  Legacy rows with
    no tracked ``ReadID`` fall back to the summed input Count.

    Junction metadata (for example ``Seq`` and ``offset``) come from the
    highest-priority dataframe that contains the junction. ``HasPEEvidence``
    records PE/DS support independently so choosing split-read sequence metadata
    does not discard the PE provenance used by downstream filters.
    """
    frames = []
    for source_priority, dataframe in enumerate(dataframes):
        if dataframe is None or dataframe.empty:
            continue

        dataframe = _upgrade_read_evidence(dataframe)
        for column in BP_PAIR_COLUMNS:
            if column not in dataframe.columns:
                if column in {'ReadID', 'Readbarcode'}:
                    dataframe[column] = [[] for _ in range(dataframe.shape[0])]
                elif column == 'HasPEEvidence':
                    dataframe[column] = False
                else:
                    dataframe[column] = None
        dataframe['_JunctionKey'] = dataframe.apply(_canonical_junction_key, axis=1)
        dataframe['_SourcePriority'] = source_priority
        frames.append(dataframe)

    if not frames:
        return pd.DataFrame(columns=BP_PAIR_COLUMNS)

    evidence = pd.concat(frames, ignore_index=True)
    merged_rows = []

    for junction_key, group in evidence.groupby('_JunctionKey', sort=False):
        # Final junction-level deduplication.  A QNAME can be observed in more
        # than one evidence source (split, PE, and DS), and the two alignments
        # of one discordant read pair also share the same stable ReadID.  Keep
        # one occurrence per ReadID while retaining the first non-empty barcode.
        unique_read_evidence = {}

        for _, row in group.iterrows():
            read_ids = _read_id_array(row['ReadID'])
            readbarcodes = _evidence_list(row['Readbarcode'])
            for index, read_id in enumerate(read_ids):
                barcode = readbarcodes[index] if index < len(readbarcodes) else None
                barcode = str(barcode) if barcode is not None and str(barcode) not in {'', 'None', 'nan'} else None
                read_id = int(read_id)
                if read_id not in unique_read_evidence:
                    unique_read_evidence[read_id] = barcode
                elif unique_read_evidence[read_id] is None and barcode is not None:
                    unique_read_evidence[read_id] = barcode

        counts = pd.to_numeric(group['Count'], errors='coerce').fillna(0)
        preferred = group[group['_SourcePriority'] == group['_SourcePriority'].min()]
        preferred_counts = pd.to_numeric(preferred['Count'], errors='coerce').fillna(0)
        representative = preferred.loc[preferred_counts.idxmax()].copy()
        end1, end2 = junction_key
        read_ids = np.fromiter(
            unique_read_evidence.keys(),
            dtype=np.uint64,
            count=len(unique_read_evidence),
        )

        representative['Chrom1'], representative['Coord1'], representative['Clip1'] = end1
        representative['Chrom2'], representative['Coord2'], representative['Clip2'] = end2
        representative['ReadID'] = read_ids
        readbarcodes = list(unique_read_evidence.values())
        if not any(barcode is not None for barcode in readbarcodes):
            readbarcodes = []

        # BARCODE: keep WGS rows sparse after ReadID deduplication.
        representative['Readbarcode'] = readbarcodes
        raw_count = int(counts.sum())
        representative['Count'] = int(read_ids.size) if read_ids.size else raw_count
        representative['HasPEEvidence'] = bool(
            group['Seq'].isin(['PE_Support', 'DS']).any()
            or group['HasPEEvidence'].map(_evidence_flag).any()
        )
        merged_rows.append(representative[BP_PAIR_COLUMNS].tolist())

    result = pd.DataFrame(merged_rows, columns=BP_PAIR_COLUMNS)
    return result.sort_values('Count', ascending=False).reset_index(drop=True)


def _build_cn_boundary_components(cn_amp, min_cn=4.5, max_gap=30000, min_boundary_clip=5):
    """Build continuous clip-derived CN components and supported boundaries."""
    work = cn_amp.loc[pd.to_numeric(cn_amp['CN'], errors='coerce').ge(min_cn)].copy()
    if 'SegSource' in work.columns:
        work = work.loc[work['SegSource'].fillna('clip').eq('clip')]

    components, boundary_rows = [], []
    component_id = 0
    for chrom, group in work.groupby('Chrom', sort=False):
        group = group.sort_values(['Start', 'End'])
        current_start = current_end = None
        current_rows = []

        def finish_component():
            nonlocal component_id, current_start, current_end, current_rows
            if not current_rows:
                return
            components.append({'Chrom': str(chrom), 'Start': int(current_start), 'End': int(current_end), 'Component': component_id})
            for segment in current_rows:
                for coord_col, clip_col, clip in (('Start', 'ClipLeft', 'L'), ('End', 'ClipRight', 'R')):
                    value = pd.to_numeric(segment.get(clip_col, 0), errors='coerce')
                    support = 0 if pd.isna(value) else int(value)
                    if support >= min_boundary_clip:
                        boundary_rows.append({'Chrom': str(chrom), 'Coord': int(segment[coord_col]), 'Clip': clip, 'ClipDepth': support, 'Component': component_id})
            component_id += 1

        for _, segment in group.iterrows():
            start, end = int(segment['Start']), int(segment['End'])
            if current_start is None or start - current_end <= max_gap:
                if current_start is None:
                    current_start, current_end = start, end
                else:
                    current_end = max(current_end, end)
                current_rows.append(segment)
            else:
                finish_component()
                current_start, current_end, current_rows = start, end, [segment]
        finish_component()

    boundary_columns = ['Chrom', 'Coord', 'Clip', 'ClipDepth', 'Component']
    component_columns = ['Chrom', 'Start', 'End', 'Component']
    boundaries = pd.DataFrame(boundary_rows, columns=boundary_columns).drop_duplicates(['Chrom', 'Coord', 'Clip'], keep='first')
    components = pd.DataFrame(components, columns=component_columns)
    return boundaries, components


def snap_ds_to_cn_boundary(bp_duo, bp_cand_stats, cn_amp, bin_norm, bp_all, min_cn=4.5, min_reads=2, min_boundary_clip=5, max_snap_distance=30000, max_cn_gap=30000, min_continuous_cn_frac=1.0):
    """Snap one unsupported intra-chromosomal DS endpoint to a nearby supported CN boundary.

    The function updates bp_duo coordinates and bp_cand_stats directly, then
    returns the adjusted bp_duo.
    """
    audit_columns = ['RowIndex', 'Chrom', 'Count', 'AnchorCoord', 'AnchorClip', 'OldCoord', 'NewCoord', 'OtherClip', 'Distance', 'Component', 'CNBinCount', 'CNMin', 'CNAmpFraction']
    if bp_duo.empty or cn_amp.empty or bin_norm is None or bin_norm.empty:
        return bp_duo

    counts = pd.to_numeric(bp_duo['Count'], errors='coerce').fillna(0)
    eligible_mask = bp_duo['Seq'].eq('DS') & counts.ge(min_reads) & bp_duo['Chrom1'].astype(str).eq(bp_duo['Chrom2'].astype(str))
    eligible = bp_duo.loc[eligible_mask]
    if eligible.empty:
        return bp_duo

    boundaries, components = _build_cn_boundary_components(cn_amp, min_cn, max_cn_gap, min_boundary_clip)
    if boundaries.empty or components.empty:
        return bp_duo

    boundary_lookup = {(str(row.Chrom), int(row.Coord), str(row.Clip)): row for row in boundaries.itertuples(index=False)}
    component_lookup = {int(row.Component): row for row in components.itertuples(index=False)}
    by_chrom_clip_component = {(str(chrom), str(clip), int(component)): group.sort_values('Coord').reset_index(drop=True) for (chrom, clip, component), group in boundaries.groupby(['Chrom', 'Clip', 'Component'], sort=False)}

    structural_proposals, endpoints_by_row = {}, {}
    for row in eligible.itertuples():
        endpoints = [(str(row.Chrom1), int(row.Coord1), str(row.Clip1)), (str(row.Chrom2), int(row.Coord2), str(row.Clip2))]
        endpoints_by_row[row.Index] = endpoints
        proposals = []
        for anchor_index, other_index in ((0, 1), (1, 0)):
            anchor, other = endpoints[anchor_index], endpoints[other_index]
            anchor_boundary = boundary_lookup.get(anchor)
            if anchor_boundary is None:
                continue
            component_id = int(anchor_boundary.Component)
            candidate_df = by_chrom_clip_component.get((other[0], other[2], component_id))
            if candidate_df is None or candidate_df.empty:
                continue
            distances = (candidate_df['Coord'] - other[1]).abs()
            nearest_distance = int(distances.min())
            nearest = candidate_df.loc[distances.eq(nearest_distance)]
            if nearest_distance > max_snap_distance or len(nearest) != 1:
                continue
            component = component_lookup[component_id]
            if not int(component.Start) - max_snap_distance <= other[1] <= int(component.End) + max_snap_distance:
                continue
            proposals.append((nearest_distance, anchor_index, other_index, anchor_boundary, nearest.iloc[0]))
        if proposals:
            structural_proposals[row.Index] = proposals

    if not structural_proposals:
        return bp_duo

    needed_by_clip = {'L': set(), 'R': set()}
    for row_index, proposals in structural_proposals.items():
        for proposal in proposals:
            needed_by_clip[endpoints_by_row[row_index][proposal[2]][2]].add(endpoints_by_row[row_index][proposal[2]])

    raw_clip_set = set()
    for clip, side in (('L', 'left'), ('R', 'right')):
        needed = needed_by_clip[clip]
        if not needed:
            continue
        needed_coords = {coord for _, coord, _ in needed}
        raw_clip = bp_all[side]
        matches = raw_clip.loc[raw_clip['Coord'].isin(needed_coords), ['Chrom', 'Coord', 'Clip']]
        raw_clip_set.update((str(chrom), int(coord), str(direction)) for chrom, coord, direction in matches.itertuples(index=False, name=None) if str(direction) == clip)

    bin_work = bin_norm.loc[:, ['Chrom', 'Coord', 'CN']].copy()
    bin_work['Chrom'] = bin_work['Chrom'].astype(str)
    bin_work['Coord'] = pd.to_numeric(bin_work['Coord'], errors='coerce')
    bin_work['CN'] = pd.to_numeric(bin_work['CN'], errors='coerce')
    cn_by_chrom = {}
    for chrom, group in bin_work.dropna(subset=['Coord']).groupby('Chrom', sort=False):
        group = group.sort_values('Coord')
        cn_by_chrom[str(chrom)] = (group['Coord'].to_numpy(dtype=np.int64), group['CN'].to_numpy(dtype=float))

    audit, used_boundaries = [], []
    eligible_lookup = {row.Index: row for row in eligible.itertuples()}
    for row_index, proposals in structural_proposals.items():
        row = eligible_lookup[row_index]
        endpoints = endpoints_by_row[row_index]
        accepted = []
        for nearest_distance, anchor_index, other_index, anchor_boundary, candidate in proposals:
            other = endpoints[other_index]
            if other in raw_clip_set:
                continue
            chrom_bins = cn_by_chrom.get(other[0])
            if chrom_bins is None:
                continue
            span_start = min(int(anchor_boundary.Coord), int(candidate.Coord))
            span_end = max(int(anchor_boundary.Coord), int(candidate.Coord))
            bin_coords, bin_cn = chrom_bins
            left = np.searchsorted(bin_coords, span_start, side='left')
            right = np.searchsorted(bin_coords, span_end, side='right')
            span_cn = bin_cn[left:right]
            span_cn = span_cn[~np.isnan(span_cn)]
            if not len(span_cn):
                continue
            amp_fraction = float((span_cn >= min_cn).mean())
            if amp_fraction < min_continuous_cn_frac:
                continue
            accepted.append((nearest_distance, anchor_index, other_index, anchor_boundary, candidate, int(len(span_cn)), float(span_cn.min()), amp_fraction))

        if not accepted:
            continue
        accepted.sort(key=lambda item: item[0])
        if len(accepted) > 1 and accepted[0][0] == accepted[1][0]:
            continue

        distance, anchor_index, other_index, anchor, candidate, cn_bin_count, cn_min, cn_amp_fraction = accepted[0]
        old_coord = endpoints[other_index][1]
        bp_duo.at[row_index, 'Coord1' if other_index == 0 else 'Coord2'] = int(candidate.Coord)
        used_boundaries.append({'Chrom': str(candidate.Chrom), 'Coord': int(candidate.Coord), 'Clip': str(candidate.Clip), 'ClipDepth': int(candidate.ClipDepth)})
        used_boundaries.append({'Chrom': str(anchor.Chrom), 'Coord': int(anchor.Coord), 'Clip': str(anchor.Clip), 'ClipDepth': int(anchor.ClipDepth)})
        audit.append({'RowIndex': row_index, 'Chrom': endpoints[other_index][0], 'Count': int(row.Count), 'AnchorCoord': int(anchor.Coord), 'AnchorClip': anchor.Clip, 'OldCoord': old_coord, 'NewCoord': int(candidate.Coord), 'OtherClip': candidate.Clip, 'Distance': distance, 'Component': int(candidate.Component), 'CNBinCount': cn_bin_count, 'CNMin': cn_min, 'CNAmpFraction': cn_amp_fraction})

    if used_boundaries:
        used = pd.DataFrame(used_boundaries)[['Chrom', 'Coord', 'Clip', 'ClipDepth']].drop_duplicates(['Chrom', 'Coord', 'Clip'])
        used['CleanBP'], used['InDepth'], used['OutDepth'] = False, 0, 0
        used = used[['Chrom', 'Coord', 'Clip', 'CleanBP', 'ClipDepth', 'InDepth', 'OutDepth']]
        if not len(bp_cand_stats.columns):
            for column in used.columns:
                bp_cand_stats[column] = pd.Series(dtype=used[column].dtype)
        bp_cand_stats.reset_index(drop=True, inplace=True)
        existing = set(bp_cand_stats[['Chrom', 'Coord', 'Clip']].itertuples(index=False, name=None))
        for boundary in used.itertuples(index=False):
            key = (boundary.Chrom, boundary.Coord, boundary.Clip)
            if key not in existing:
                bp_cand_stats.loc[len(bp_cand_stats)] = {column: getattr(boundary, column, np.nan) for column in bp_cand_stats.columns}
                existing.add(key)

    audit_df = pd.DataFrame(audit, columns=audit_columns)
    print(f'DS endpoints snapped to CN boundaries: {len(audit_df)}')
    if not audit_df.empty:
        print(audit_df.to_string(index=False))
    return bp_duo
