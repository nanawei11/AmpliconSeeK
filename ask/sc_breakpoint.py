"""ALL raw-pair matrices for A/(A+B+1), independent of cell/ecDNA assignment."""
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd

try:
    from . import sc_matrix as sc
except ImportError:
    import sc_matrix as sc

PAIR_COLUMNS = ['Chrom1', 'Coord1', 'Clip1', 'Chrom2', 'Coord2', 'Clip2']
ANNOTATION_COLUMNS = ['JunctionID'] + PAIR_COLUMNS + [
    'MaxRawCount', 'RawRecordCount', 'Endpoint1AmpliconIDs',
    'Endpoint2AmpliconIDs', 'SharedAmpliconIDs', 'BothEndpointsInEcDNA']


def pair_key(row):
    return tuple(sorted((str(row[f'Chrom{i}']), int(row[f'Coord{i}']),
                         str(row[f'Clip{i}'])) for i in (1, 2)))


def annotate_pairs(raw, circular, min_count=5):
    """Exact, orientation-independent pairs; closed ASK segment boundaries.

    Both ends must belong to the union of actual segments, possibly different
    circles. Gaps between segments never count. Count threshold is per raw row.
    """
    records = defaultdict(list)
    for _, row in raw.iterrows():
        records[pair_key(row)].append(row)
    segments = defaultdict(list)
    for row in circular.itertuples(index=False):
        segments[str(row.Chrom)].append((int(row.Start), int(row.End), str(row.AmpliconID)))
    rows = []
    for key, entries in sorted(records.items()):
        maximum = max(float(row['Count']) for row in entries)
        if maximum < min_count:
            continue
        memberships = [{amp for start, end, amp in segments[chrom] if start <= pos <= end}
                       for chrom, pos, clip in key]
        rows.append([sc.junction_id(*key[0], *key[1]), *key[0], *key[1], maximum,
                     len(entries), ';'.join(sorted(memberships[0])),
                     ';'.join(sorted(memberships[1])),
                     ';'.join(sorted(memberships[0] & memberships[1])),
                     bool(memberships[0] and memberships[1])])
    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)


def normal_endpoint_counts(bam, chrom, pos, allowed, mapq=20):
    names = defaultdict(set)
    for read in bam.fetch(chrom, max(0, pos), pos + 1):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.mapping_quality < mapq or not read.has_tag('NM')
                or read.get_tag('NM') > 1 or not read.cigartuples
                or any(op in (4, 5) for op, length in read.cigartuples)
                or not read.has_tag('CB')):
            continue
        barcode = sc.clean_barcode(read.get_tag('CB'))
        if barcode and (allowed is None or barcode in allowed) and any(
                start <= pos < end for start, end in read.get_blocks()):
            names[barcode].add(read.query_name)
    return Counter({barcode: len(values) for barcode, values in names.items()})


def build_all_breakpoint_matrices(raw, circular, bam_path, outprefix=None,
                                  cell_barcodes=None, min_count=5, mapq=20):
    """Return support, normal, ALL annotation, ecDNA-union annotation.

    No positive-cell filtering. Without a whitelist, retain every raw support
    barcode and every normal-coverage CB observed at selected endpoints.
    Normal counts require an indexed BAM; missing evidence is never faked as 0.
    """
    raw = sc.clean_barcode_columns(sc.read_table_or_df(raw))
    circular = sc.read_table_or_df(circular)
    required = set(PAIR_COLUMNS + ['Count'])
    if not required.issubset(raw.columns):
        raise ValueError(f'Raw pair table missing columns: {sorted(required - set(raw.columns))}')
    barcode_col = next((c for c in ('Readbarcode', 'Readsbarcode') if c in raw.columns), None)
    if barcode_col is None:
        raise ValueError('Raw pair table requires Readbarcode or Readsbarcode')
    if bam_path is None or sc.pysam is None:
        raise ValueError('ALL normal matrix requires pysam and an indexed BAM')
    allowed = None if cell_barcodes is None else {
        b for value in cell_barcodes if (b := sc.clean_barcode(value))}
    annotation = annotate_pairs(raw, circular, min_count)
    selected = {
        pair_key(row): row['JunctionID'] for _, row in annotation.iterrows()}
    support = defaultdict(Counter)
    cells = set() if allowed is None else set(allowed)
    for _, row in raw.iterrows():
        barcodes = sc.clean_barcode_list(row[barcode_col])
        barcodes = [b for b in barcodes if allowed is None or b in allowed]
        cells.update(barcodes)
        key = pair_key(row)
        if key in selected:
            support[selected[key]].update(barcodes)
    endpoints = sorted({(chrom, pos) for key in selected for chrom, pos, clip in key})
    normal = {}
    with sc.pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if not bam.has_index():
            raise ValueError('ALL normal matrix requires an indexed BAM')
        missing = {chrom for chrom, pos in endpoints} - set(bam.references)
        if missing:
            raise ValueError(f'Raw breakpoint chromosomes absent from BAM: {sorted(missing)}')
        for chrom, pos in endpoints:
            normal[(chrom, pos)] = normal_endpoint_counts(bam, chrom, pos, allowed, mapq)
            cells.update(normal[(chrom, pos)])
    cells = sorted(cells)
    ids = annotation.JunctionID.tolist()
    a = pd.DataFrame(0, index=pd.Index(cells, name='Barcode'), columns=ids)
    b = a.copy()
    for key, jid in selected.items():
        for barcode, count in support[jid].items():
            a.at[barcode, jid] = count
        counts = normal[key[0][:2]] + normal[key[1][:2]]
        for barcode, count in counts.items():
            b.at[barcode, jid] = count
    a, b = a.reset_index(), b.reset_index()
    ecdna = annotation.loc[annotation.BothEndpointsInEcDNA.astype(bool)].copy()
    if outprefix is not None:
        Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
        for suffix, table in [('cell_by_breakpoint_support', a), ('cell_by_breakpoint_normal', b),
                              ('breakpoint_annotation_ALL', annotation),
                              ('breakpoint_annotation_ecDNA', ecdna)]:
            table.to_csv(f'{outprefix}_{suffix}.tsv', sep='\t', index=False)
    return a, b, annotation, ecdna


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', required=True, help='ASK breakpoint_pair_raw.tsv')
    parser.add_argument('--circular', required=True, help='ASK circular segment TSV')
    parser.add_argument('--derived-circular', help='Optional derived circle segments')
    parser.add_argument('--bam', required=True, help='Original indexed BAM')
    parser.add_argument('--outprefix', required=True, help='Prefix ending in _ask_sc')
    parser.add_argument('--barcodes', help='One valid-cell barcode per line; no header')
    parser.add_argument('--min-count', type=int, default=5)
    args = parser.parse_args()
    circular = sc.read_table_or_df(args.circular)
    if args.derived_circular:
        circular = pd.concat([circular, sc.read_table_or_df(args.derived_circular)], ignore_index=True)
    out = Path(args.outprefix).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    method = out.parent / '\u505a\u6cd5.md'
    config = '\n\n## ALL breakpoint matrix export: ' + out.name + '\n'
    config += f'Generator: {Path(__file__).resolve()}; dependency: {Path(sc.__file__).resolve()}\n'
    config += '\n'.join(f'{key}: {Path(value).resolve()}' for key, value in vars(args).items()
                        if key in ('raw','circular','derived_circular','bam','barcodes') and value)
    config += f'\nOutput prefix: {out}\nCount threshold: {args.min_count}; MAPQ>=20, NM<=1; status: running.\n'
    with method.open('a') as handle:
        handle.write(config)
    try:
        a,b,ann,ec = build_all_breakpoint_matrices(args.raw,circular,args.bam,out,
            cell_barcodes=Path(args.barcodes).read_text().splitlines() if args.barcodes else None,
            min_count=args.min_count)
    except Exception as error:
        with method.open('a') as handle:
            handle.write(f'Status: failed; {type(error).__name__}: {error}\n')
        raise
    with method.open('a') as handle:
        handle.write(f'Status: completed; cells={len(a)}, ALL pairs={len(ann)}, ecDNA-region pairs={len(ec)}.\n'
                     'Canonicalize exact coordinates and Clip directions regardless of endpoint order. Select a pair if any raw row meets the Count threshold; '
                     'sum barcode occurrences across all raw rows for that pair. Normal coverage requires primary alignments without soft/hard clipping and with CB; '
                     'aligned blocks must cover the exact endpoint. Deduplicate QNAME per cell per endpoint, then sum both endpoints. '
                     'Use the union of closed ecDNA segment intervals; both endpoints must overlap this union and may belong to different circles. '
                     'ALL is not restricted to ecDNA-positive cells. No independent real-BAM spot check was performed.\n')
        for suffix in ('cell_by_breakpoint_support','cell_by_breakpoint_normal',
                       'breakpoint_annotation_ALL','breakpoint_annotation_ecDNA'):
            handle.write(f'Output: {out}_{suffix}.tsv\n')
