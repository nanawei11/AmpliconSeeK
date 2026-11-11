################################################################################
# breakpoint pair detection and filtering
################################################################################
#------------------------------------------------------------------------------#
import pysam
import re
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
from grange import GRange
import misc


#------------------------------------------------------------------------------#
# filter breakpoint pairs
#------------------------------------------------------------------------------#
# 2026-0831-adj
_SEQ_PATTERN = re.compile(r'([A-Z]+)')
_EMPTY_INTERVALS = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))

# columns the filters below actually read; everything else in bp_duo is
# carried through untouched and only materialized at the end.
_WORK_COLS = ('Chrom1', 'Coord1', 'Clip1', 'Chrom2', 'Coord2', 'Clip2',
              'Count', 'offset', 'Seq', 'HasPEEvidence')


def get_match_seq(seq):
    match = _SEQ_PATTERN.search(str(seq))
    return match.group(0) if match else ''


def _merge_intervals_by_chrom(cn_amp, cn_amp_error):
    """Pad every cn_amp interval by +/- cn_amp_error, then collapse each
    chromosome's intervals into sorted, non-overlapping runs.

    Padding can make neighboring segments overlap even when the raw
    segments don't, so merging first is what lets a single binary search
    answer membership correctly.
    """
    merged = {}
    for chrom, group in cn_amp.groupby('Chrom'):
        starts = group['Start'].to_numpy() - cn_amp_error
        ends = group['End'].to_numpy() + cn_amp_error
        order = np.argsort(starts, kind='mergesort')
        starts, ends = starts[order], ends[order]

        merged_starts = [starts[0]]
        merged_ends = [ends[0]]
        for s, e in zip(starts[1:], ends[1:]):
            if s <= merged_ends[-1]:
                if e > merged_ends[-1]:
                    merged_ends[-1] = e
            else:
                merged_starts.append(s)
                merged_ends.append(e)

        merged[chrom] = (np.asarray(merged_starts), np.asarray(merged_ends))
    return merged


def _in_intervals(chrom_arr, coord_arr, merged_intervals):
    """Vectorized membership: is each (chrom, coord) inside one of that
    chromosome's merged, padded cn_amp intervals?

    Returns a positionally-aligned bool ndarray. No temporary DataFrame --
    chromosomes are grouped with factorize + boolean masks.
    """
    chrom_arr = np.asarray(chrom_arr)
    coord_arr = np.asarray(coord_arr)
    result = np.zeros(chrom_arr.size, dtype=bool)
    if chrom_arr.size == 0:
        return result

    codes, uniques = pd.factorize(chrom_arr)
    for i, chrom in enumerate(uniques):
        starts, ends = merged_intervals.get(chrom, _EMPTY_INTERVALS)
        if starts.size == 0:
            continue
        sel = np.flatnonzero(codes == i)
        coords = coord_arr[sel]
        idx = np.searchsorted(starts, coords, side='right') - 1
        valid = idx >= 0
        hit = np.zeros(coords.size, dtype=bool)
        hit[valid] = coords[valid] <= ends[idx[valid]]
        result[sel] = hit
    return result


def bp_pair_filter(bp_duo, cn_amp=None, bp_cand_stats=None,
                   short_distance=200, min_junc_cnt=5,
                   max_offset_overhang=1000, min_match_len=5,
                   max_offset_insertion=20, max_srp=0.7,
                   cn_amp_error=10 * 1000):

    print(f'initial bp pair: {bp_duo.shape}')

    n_col = bp_duo.shape[1]
    work = bp_duo[[c for c in _WORK_COLS if c in bp_duo.columns]]
    seq_all = bp_duo['Seq'].to_numpy()

    # 1. Filter by Count (positional survivors only -- no wide-table copy)
    pos = np.flatnonzero(bp_duo['Count'].to_numpy() >= min_junc_cnt)
    print(f'after Count >= {min_junc_cnt}: ({pos.size}, {n_col})')

    stage = work.iloc[pos]

    # 2. Filter short same-chromosome distances and offsets
    short_fail = (
        stage['Chrom1'].eq(stage['Chrom2'])
        & stage['Coord2'].sub(stage['Coord1']).abs().lt(short_distance)
    )

    offset_fail = (
        stage['offset'].gt(max_offset_insertion)
        | stage['offset'].lt(-max_offset_overhang)
    )

    basic_fail = short_fail | offset_fail

    print(f'remove based on short_distance: {short_fail.sum()}')
    print(f'remove additionally based on offset: {(offset_fail & ~short_fail).sum()}')

    # 3. Identify DS and PE evidence
    is_ds = stage['Seq'].eq('DS')
    pe_evidence = stage['Seq'].eq('PE_Support')

    if 'HasPEEvidence' in stage.columns:
        pe_flag = stage['HasPEEvidence'].map(
            lambda value:
                value.strip().lower() in {'1', 'true', 't', 'yes', 'y'}
                if isinstance(value, str)
                else bool(value) if pd.notna(value)
                else False
        )
        pe_evidence = pe_evidence | pe_flag

    skip_sequence_filter = is_ds | pe_evidence
    sequence_rows = ~skip_sequence_filter
    sequence_fail = pd.Series(False, index=stage.index)

    match_seq = stage.loc[sequence_rows, 'Seq'].map(get_match_seq)
    match_len = match_seq.str.len()
    srp = match_seq.map(
        lambda seq: simple_repeat_proportion(seq) if len(seq) > 0 else 0)

    sequence_fail.loc[sequence_rows] = (srp.gt(max_srp) | match_len.lt(min_match_len))

    print(f'remove based on sequence: {sequence_fail.sum()}')

    # Apply basic and sequence filters
    keep = ~basic_fail & (skip_sequence_filter | ~sequence_fail)
    pos = pos[keep.to_numpy()]

    print(f'after basic filter: ({pos.size}, {n_col})')
    print(f"remaining DS: {(seq_all[pos] == 'DS').sum()}")

    # 5. Both endpoints must fall within a cn_amp segment, allowing cn_amp_error padding
    if cn_amp is not None:

        merged_intervals = _merge_intervals_by_chrom(cn_amp, cn_amp_error)
        stage = work.iloc[pos]

        end1_in_amp = _in_intervals(
            stage['Chrom1'].to_numpy(), stage['Coord1'].to_numpy(), merged_intervals)

        # End2 only matters where End1 already passed (filter is End1 AND End2)
        cn_keep = np.zeros(pos.size, dtype=bool)
        sel = np.flatnonzero(end1_in_amp)
        if sel.size:
            cn_keep[sel] = _in_intervals(
                stage['Chrom2'].to_numpy()[sel],
                stage['Coord2'].to_numpy()[sel],
                merged_intervals)

        print(f'remove based on cn_amp: {(~cn_keep).sum()}')
        pos = pos[cn_keep]
        print(f'after cn_amp filter: ({pos.size}, {n_col})')

    # 6. DS requires at least one endpoint in the breakpoint table
    if bp_cand_stats is not None:

        stage = work.iloc[pos]
        is_ds = (seq_all[pos] == 'DS')

        # Retain non-DS rows; query the breakpoint table only for DS rows
        breakpoint_keep = np.ones(pos.size, dtype=bool)
        ds_sel = np.flatnonzero(is_ds)
        if ds_sel.size:
            bp_index = pd.MultiIndex.from_frame(
                bp_cand_stats[['Chrom', 'Coord', 'Clip']].drop_duplicates())
            ds_rows = stage.iloc[ds_sel]

            end1_in_bp = pd.MultiIndex.from_arrays(
                [ds_rows['Chrom1'], ds_rows['Coord1'], ds_rows['Clip1']]).isin(bp_index)
            end2_in_bp = pd.MultiIndex.from_arrays(
                [ds_rows['Chrom2'], ds_rows['Coord2'], ds_rows['Clip2']]).isin(bp_index)

            breakpoint_keep[ds_sel] = end1_in_bp | end2_in_bp

        print(
            'remove DS with neither endpoint in breakpoint: '
            f'{(~breakpoint_keep).sum()}'
        )

        pos = pos[breakpoint_keep]

    # the only wide-table materialization, on the smallest row set
    bp_pair = bp_duo.iloc[pos].copy()

    print(f'final bp pair: {bp_pair.shape}')
    print(f"final DS: {bp_pair['Seq'].eq('DS').sum()}")

    return bp_pair.reset_index(drop=True)

# def bp_pair_filter(bp_duo, short_distance = 200, min_junc_cnt = 5,
#                    max_offset_overhang = 1000, min_match_len = 5,
#                    max_offset_insertion = 20, max_srp = 0.7):
#     """
#     filter breakpoint pairs

#     max_srp : max simple repeat proportion
#     """
#     print(f'inital bp pair: {bp_duo.shape}')
#     bp_pair = bp_duo[bp_duo['Count'] >= min_junc_cnt]
#     print(f'inital bp pair: {bp_pair.shape}')

#     basic_fail = (bp_pair['Chrom1'] == bp_pair['Chrom2']) \
#         & (abs(bp_pair['Coord2'] - bp_pair['Coord1']) < short_distance)
#     print(f'remove based on short_distance: {sum(basic_fail)}')

#     # also control for max insertion size (+) and max overhang size (-)
#     basic_fail = basic_fail | (bp_pair['offset'] > max_offset_insertion) \
#                             | (bp_pair['offset'] < -max_offset_overhang)
#     print(f'remove based on max insertion size: {sum(basic_fail)}')

#     # DS rows have no junction matching sequence, so do not apply the
#     # repeat/match-length filter to them. They still have to pass Count,
#     # short-distance and offset filters above.
#     is_ds = bp_pair['Seq'].eq('DS')
#     sequence_fail = pd.Series(False, index=bp_pair.index)
#     sequence_rows = ~is_ds
#     match_seq = bp_pair.loc[sequence_rows, 'Seq'].map(
#         lambda seq: re.search('([A-Z]+)', seq).group(0))
#     match_len = match_seq.map(len)
#     srp = match_seq.map(simple_repeat_proportion)
#     sequence_fail.loc[sequence_rows] = \
#         (srp > max_srp) | (match_len < min_match_len)

#     print(f'remove based on simple repeat: {sum(sequence_fail)}')
#     # Keep the existing sequence-filter exemption for other PE-supported rows.
#     pe_evidence = bp_pair['Seq'].eq('PE_Support')
#     if 'HasPEEvidence' in bp_pair.columns:
#         flag = bp_pair['HasPEEvidence'].map(
#             lambda value: value.strip().lower() in {'1', 'true', 't', 'yes', 'y'}
#             if isinstance(value, str)
#             else bool(value) if pd.notna(value) else False
#         )
#         pe_evidence = pe_evidence | flag
#     skip_sequence_filter = is_ds | pe_evidence
#     keep = ~basic_fail & (skip_sequence_filter | ~sequence_fail)
#     print(f'after filter: {bp_pair.loc[keep].shape}')
#     return bp_pair.loc[keep].reset_index(drop=True)

def estimaed_thred(bp_duo, i='current'):
    '''
    Estimated thred used for bp pair filtered
    '''
    if (bp_duo.shape[0] > 10):
        try:
            data_ = np.log(bp_duo['Count'])
        except:
            data_ = np.log(bp_duo['Count'].astype('int'))
        counts_, bins_ = np.histogram(data_, bins=1000)
        nozero = counts_[counts_!=0]
        # diff_ratio = (nozero[:-1] - nozero[1:])/nozero[1:]
        # print(diff_ratio[:10])
        # all_dat[i] = {}
        # all_dat[i]['bpduo_count'] = bp_duo['Count']
        # all_dat[i]['data_raw'] = data_
        # all_dat[i]['counts_'] = counts_
        # all_dat[i]['nozero'] = nozero
        # all_dat[i]['bin_values'] = bins_
        # all_dat[i]['diff_ratio'] = diff_ratio
        cusum_ratio = (nozero[:-1] / np.cumsum(nozero)[1:])

        if nozero.shape[0] <= 5:
            index = -1
        else:
            index = np.argmax(cusum_ratio <= 0.01)

        if (index >= 9) and (bp_duo.shape[0] <= 10000) and (index > 0):
            if bp_duo.shape[0] <= 2000:
                index = np.min([np.argmax(cusum_ratio <= 0.05), 4])
            elif bp_duo.shape[0] <= 4000:
                index = np.min([np.argmax(cusum_ratio <= 0.05), 6])
            else:
                index = np.min([np.argmax(cusum_ratio <= 0.05), 10])
                
        idx = np.argmax(counts_ == nozero[index+1])
        thred = np.exp(bins_[idx:idx+2].mean())

        filter_ = bp_duo['Count'][ bp_duo['Count'] >= round(thred,0) ]

        if filter_.shape[0] > 30000 :
            thred = np.exp(data_.sort_values(ascending=False)[1200])    

        if (bp_duo.shape[0] > 10000) and filter_.shape[0] > 1000:
            thred = np.max([thred, 20])

        if (bp_duo.shape[0] <= 5000) and filter_.shape[0] < 300:
            thred = np.min([thred, 5])
            
        if (bp_duo.shape[0] > 1000) and filter_.shape[0] >= 300:
            thred = np.max([thred, 3])

        filter_ = bp_duo['Count'][ bp_duo['Count'] >= round(thred,0) ]
        print(f'sample:{i}, index: {index}, thred: {thred}, data.shape: {bp_duo.shape[0]}, filter shape: {filter_.shape[0]}')
        return index, thred, bp_duo.shape[0], filter_.shape[0]
    
#------------------------------------------------------------------------------#
def simple_repeat_proportion(s, n_iter = 4):
    """search simple repeat in a string
    return the proportion
    """
    match_ = [i for i in re.compile(r"(.+?)\1+").finditer(s)]
    return sum([len(i.group()) for i in match_
                if len(i.group()) >= n_iter])/len(s)

#------------------------------------------------------------------------------#
# get breakpoint from bp pairs and cn segments
#------------------------------------------------------------------------------#
def prepare_segment_for_graph(cn_amp):
    """
    merge normal clip segments, but keep PE alternative segments independently
    """

    if 'SegSource' not in cn_amp.columns:
        return pd.DataFrame( misc.merge_bed(cn_amp), columns=['Chrom', 'Start', 'End'])

    #  raw Clip segment keep merge
    clip_seg = cn_amp[cn_amp['SegSource'] == 'clip']
    if clip_seg.empty:
        clip_merged = pd.DataFrame(columns=['Chrom', 'Start', 'End'])
    else:
        clip_merged = pd.DataFrame( misc.merge_bed(clip_seg), columns=['Chrom', 'Start', 'End'] )

    # PE segment do not merge, as inner edge
    pe_seg = cn_amp.loc[ cn_amp['SegSource'].isin(['PE_boundary', 'PE_seg']),
        ['Chrom', 'Start', 'End']
    ].copy()

    seg = pd.concat([clip_merged, pe_seg])
    seg = seg.drop_duplicates(['Chrom', 'Start', 'End'])

    return seg.sort_values(['Chrom', 'Start', 'End']).reset_index(drop=True)


def snap_pe_to_segment_boundary(bp_duo, cn_amp, tolerance=15):
    """Snap PE breakpoints near amplified-segment boundaries to boundary coordinates."""
    out = bp_duo.copy().reset_index(drop=True)
    cn_amp_graph = prepare_segment_for_graph(cn_amp)
    boundaries = {'L': {}, 'R': {}}

    for row in cn_amp_graph.itertuples(index=False):
        chrom, start, end = str(row.Chrom), int(row.Start), int(row.End)
        boundaries['L'].setdefault(chrom, []).append(start)
        boundaries['R'].setdefault(chrom, []).append(end)

    def nearest_boundary(chrom, coord, clip):
        candidates = boundaries.get(clip, {}).get(chrom, [])
        nearby = [(abs(int(coord) - candidate), candidate) for candidate in candidates if abs(int(coord) - candidate) <= tolerance]
        if not nearby:
            return None

        nearby.sort()
        if len(nearby) > 1 and nearby[0][0] == nearby[1][0]:
            return None
        return nearby[0][1]

    for index, row in out.iterrows():
        if row.get('Seq') != 'PE_Support':
            continue

        coord1 = nearest_boundary(row['Chrom1'], row['Coord1'], row['Clip1'])
        coord2 = nearest_boundary(row['Chrom2'], row['Coord2'], row['Clip2'])
        if coord1 is None or coord2 is None:
            continue

        out.at[index, 'Coord1'] = coord1
        out.at[index, 'Coord2'] = coord2

    return out


def bp_refine(bp_pair, bp_cand_stats, cn_amp, min_cn = 5):
    """
    Get breakpoints using amplified segments whose rounded CN >= min_cn.
    Round a local copy so callers retain their original CN estimates.
    """
    cn_amp = cn_amp.copy()
    cn_amp['CN'] = pd.to_numeric(cn_amp['CN'], errors='raise').round()
    cn_amp = cn_amp[cn_amp['CN'] >= min_cn]
    # add sequence depth to stats
    col1 = ['Chrom1', 'Coord1', 'Clip1']
    col2 = ['Chrom2', 'Coord2', 'Clip2']
    newcol = ['Chrom', 'Coord', 'Clip']
    bp_df = pd.concat([bp_pair[col1].set_axis(newcol, axis=1), \
                    bp_pair[col2].set_axis(newcol, axis=1)])
    bp_df = bp_df.drop_duplicates()

    # add clip depth to stats
    bp_stats = pd.merge(bp_df, bp_cand_stats,
                        on=['Chrom', 'Coord', 'Clip'], how='left')

    # breakpoints from amplified segments
    ## merge adjacent segments into one
    # cn_amp_merged = misc.merge_bed(cn_amp)
    graph_seg = prepare_segment_for_graph(cn_amp)
    # Only original clip-derived segments may provide ClipDepth provenance.
    # PE_boundary/PE_seg rows can contain nearby clip annotations, but their
    # boundaries were introduced by PE evidence and must not be reclassified
    # as clip breakpoints for get_segment().
    cn_amp_clip = dict()
    cn_amp_disc = dict()
    for i in cn_amp.itertuples():
        if str(getattr(i, 'SegSource', 'clip')) == 'clip':
            cn_amp_clip[(i.Chrom, i.Start, 'L')] = i.ClipLeft
            cn_amp_clip[(i.Chrom, i.End, 'R')] = i.ClipRight
        cn_amp_disc[(i.Chrom, i.Start, 'L')] = getattr(i, 'DiscLeft', 0) or 0
        cn_amp_disc[(i.Chrom, i.End, 'R')] = getattr(i, 'DiscRight', 0) or 0

    # A PE edge can land on exactly the same coordinate as an original clip
    # segment boundary.  Restore that boundary's clip annotation after the
    # left merge, but never substitute PE support for ClipDepth.
    missing_clip = bp_stats['ClipDepth'].isna()
    if missing_clip.any():
        bp_stats.loc[missing_clip, 'ClipDepth'] = bp_stats.loc[missing_clip].apply(
            lambda row: cn_amp_clip.get(
                (row['Chrom'], row['Coord'], row['Clip']), np.nan),
            axis=1)

    # A breakpoint may enter bplist through bp_pair support, such as discordant
    # PE_Support, without matching clip evidence in bp_cand_stats. The merge then
    # leaves CleanBP as NaN and the later filter removes the breakpoint. Keep it
    # as a clean breakpoint when its cn_amp boundary has DiscLeft/DiscRight
    # support, because missing clip evidence alone does not make it unreliable.
    need_fix = bp_stats['CleanBP'].isna()
    if need_fix.any():
        bp_stats.loc[need_fix, 'CleanBP'] = bp_stats.loc[need_fix].apply(
            lambda row: cn_amp_disc.get((row['Chrom'], row['Coord'], row['Clip']), 0) > 0,
            axis=1)

    ## make breakpoints dataframe
    bplist = set(i[1:4] for i in bp_stats.itertuples())
    op = []
    for row in graph_seg.itertuples(index=False):
        L = (row.Chrom, row.Start, 'L')
        R = (row.Chrom, row.End, 'R')
        if L not in bplist:
            op.append([row.Chrom, row.Start, 'L', True,
                       cn_amp_clip.get(L, np.nan), 0, 0])
        if R not in bplist:
            op.append([row.Chrom, row.End, 'R', True,
                       cn_amp_clip.get(R, np.nan), 0, 0])
    cn_seg_df = pd.DataFrame(op, columns=bp_stats.columns)
    # cn_seg_df

    # merge and output
    df = pd.concat([bp_stats, cn_seg_df])
    df = df[pd.notnull(df.CleanBP)]
    df = df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
    return df.reset_index(drop=True), graph_seg

#------------------------------------------------------------------------------#
# infer segments from breakpoint pairs
#------------------------------------------------------------------------------#
def get_segment_on_cnamp(cn_amp_merged, bp_pair_m, knn = 3):
    '''
    infer segments from bp pair table on cn ampilified

    usage:
        get_segment_on_cnamp(cn_amp, bp_pair)

    input:
    cn_amp_merged: amplified intervals in either the legacy list format or a
        dataframe containing Chrom/Start/End columns
    bp_pair_m: pandas dataframe of bp pair
    '''
    # FIX(20260802, DataFrame adapter): ``bp_refine`` now returns
    # ``prepare_segment_for_graph`` as a dataframe,
    # whereas the original implementation returned ``misc.merge_bed`` as a
    # list.  Iterating a dataframe directly yields column names, so normalize
    # both representations to (Chrom, Start, End) tuples here.
    if isinstance(cn_amp_merged, pd.DataFrame):
        required = ['Chrom', 'Start', 'End']
        missing = [col for col in required if col not in cn_amp_merged.columns]
        if missing:
            raise ValueError(
                'cn_amp_merged dataframe is missing columns: ' +
                ', '.join(missing))
        cn_amp_regions = list(cn_amp_merged[required].itertuples(
            index = False, name = None))
    else:
        cn_amp_regions = list(cn_amp_merged)

    # record the breakpoints that have been traversed
    rm_bp_fine = []
    seg_subcn = []
    def get_subseg_R(bp_, bp_pair_df_r, rm_bp_fine, seg_subcn, knn = knn):
        # find the end point for left breakpoint
        if not bp_pair_df_r.empty:
            try:
                dst = bp_[1] - bp_pair_df_r['Coord1']
                cand = bp_pair_df_r.loc[dst < 0].reset_index(drop=True)
                ind = np.argsort(cand['Coord2'])
                # if (bp_[1] == 55904441) or (bp_[1] == 56049369):
                #     print(f'  ******* bp is L: {bp_}  from  {cand.iloc[ind[:knn],:]}*********')
                if ind.shape[0] == 1:
                    seg_subcn.append([bp_[0]] + [bp_[1]] + [cand['Coord1'][ind[0]]])
                    rm_bp_fine.append(bp_)
                else:
                    cand = cand.iloc[ind[:knn],:].reset_index(drop=True)
                    ind = np.argmax(cand['Count'])
                    diff_first_max_ratio = (cand['Count'][ind] -  cand['Count'][0]) / cand['Count'][ind]
                    if (diff_first_max_ratio <= 0.5 ) & (diff_first_max_ratio!=0):   
                        seg_subcn.append([bp_[0]] + [bp_[1]] + [cand['Coord1'][0]])
                    seg_subcn.append([bp_[0]] + [bp_[1]] + [cand['Coord1'][ind]])
                    rm_bp_fine.append(bp_)

            except:
                1
                # print(bp_)
        return  seg_subcn, rm_bp_fine
        
    def get_subseg_L(bp_, bp_pair_df_l, rm_bp_fine, seg_subcn, knn = knn):
        # find the start point for right breakpoint
        if not bp_pair_df_l.empty:
            try:
                dst = bp_pair_df_l['Coord2'] - bp_[1]
                cand = bp_pair_df_l.loc[dst < 0].reset_index(drop=True)
                ind = np.argsort(-cand['Coord1'])
                # if (bp_[1] == 55904441) or (bp_[1] == 56049369):
                #     print(f'  ******* bp is R: {bp_}  from  {cand.iloc[ind[:knn],:]}*********')
                if ind.shape[0] == 1:
                    seg_subcn.append([bp_[0]] + [cand['Coord2'][ind[0]]] + [bp_[1]])
                    rm_bp_fine.append(bp_)
                else:
                    cand = cand.iloc[ind[:knn],:].reset_index(drop=True)
                    ind = np.argmax(cand['Count'])
                    diff_first_max_ratio = (cand['Count'][ind] -  cand['Count'][0]) / cand['Count'][ind]
                    if (diff_first_max_ratio <= 0.5 ) & (diff_first_max_ratio!=0):   
                        seg_subcn.append([bp_[0]] + [cand['Coord2'][0]] + [bp_[1]])
                    seg_subcn.append([bp_[0]] + [cand['Coord2'][ind]] + [bp_[1]])
                    rm_bp_fine.append(bp_)
                # if ind.shape[0] >  2:
                #     seg_subcn.append([bp_[0]] + [bp_pair_df_l.loc[dst < 0, 'Coord2'].iloc[ind[1]]] + [bp_[1]])    
            except:
                # print(bp_)
                1
        return  seg_subcn, rm_bp_fine  

    for i in cn_amp_regions:
        if len(i) < 3:
            raise ValueError(
                'each cn_amp_merged interval must contain Chrom, Start, End')
        i = (str(i[0]), int(i[1]), int(i[2]))
        bp_pair_chr = bp_pair_m[bp_pair_m['Chrom1'] == i[0]]

        # get the candidate bp pair on segment
        ind = (bp_pair_chr['Clip1'] == 'R') & (bp_pair_chr['Clip2'] == 'L') & \
            (bp_pair_chr['Coord1'] >= i[1]) & (bp_pair_chr['Coord1']  <= i[2])
        bp_pair_df_r = bp_pair_chr[ind]

        ind = (bp_pair_chr['Clip1'] == 'R') & (bp_pair_chr['Clip2'] == 'L') & \
            (bp_pair_chr['Coord2'] >= i[1]) & (bp_pair_chr['Coord2']  <= i[2])
        bp_pair_df_l = bp_pair_chr[ind]

        # get the bp_pair on segment
        bp_pair_df  = bp_pair_chr[(bp_pair_chr['Coord1'] >= i[1]) & (bp_pair_chr['Coord1'] <= i[2]) | 
                                (bp_pair_chr['Coord2'] >= i[1]) & (bp_pair_chr['Coord2'] <= i[2])]

        for j in bp_pair_df.itertuples():
            bp_ = j[1:4]
            if bp_ not in rm_bp_fine:
                if (bp_[1] >= i[1]) and (bp_[1] <= i[2]):  # bp_ in cn_amp
                    if bp_[2] == 'L':
                        seg_subcn, rm_bp_fine  = get_subseg_R(bp_, bp_pair_df_r, rm_bp_fine, seg_subcn)
                    else:
                        seg_subcn, rm_bp_fine = get_subseg_L(bp_, bp_pair_df_l, rm_bp_fine, seg_subcn)
                
            bp_ = j[4:7]
            if bp_ not in rm_bp_fine:
                if (bp_[1] >= i[1]) & (bp_[1] <= i[2]):
                    if bp_[2] == 'L':
                        seg_subcn, rm_bp_fine = get_subseg_R(bp_, bp_pair_df_r, rm_bp_fine, seg_subcn)
                    else:
                        seg_subcn, rm_bp_fine = get_subseg_L(bp_, bp_pair_df_l, rm_bp_fine, seg_subcn)
    return rm_bp_fine, seg_subcn


def get_segment(bp_fine, bp_pair, cn_amp_merged, restrict = False, \
        min_segment_size = 100, max_segment_size = 50000000, subseg = False,
        rm_amp_df = False, knn = 3, min_component_fraction = 0.5,
        component_boundary_tolerance = 20000):
    """
    Infer genomic segments from refined breakpoint boundaries.

    Newly inferred segments are searched only within one continuous
    CN-amplified component, allowing ``component_boundary_tolerance`` bases at
    its boundaries. Inferred records must overlap merged amplified regions for
    at least ``min_component_fraction`` of their length. The boundary tolerance
    protects true segments whose breakpoint extends just beyond the CN call,
    while the overlap threshold rejects long bridges across normal regions.
    Sorted NumPy arrays and binary search replace repeated whole-chromosome
    dataframe scans, keeping working memory linear in the number of
    breakpoints.
    """
    if knn < 1:
        raise ValueError('knn must be >= 1')
    if not 0 <= min_component_fraction <= 1:
        raise ValueError('min_component_fraction must be between 0 and 1')
    if component_boundary_tolerance < 0:
        raise ValueError('component_boundary_tolerance must be >= 0')

    colnames = ['Chrom', 'Start', 'End']

    # Normalize the CN intervals once and reuse this compact table throughout.
    if isinstance(cn_amp_merged, pd.DataFrame):
        missing = [x for x in colnames if x not in cn_amp_merged.columns]
        if missing:
            raise ValueError(
                'cn_amp_merged is missing columns: ' + ', '.join(missing))
        amp_regions = cn_amp_merged.loc[:, colnames].copy()
    else:
        amp_regions = pd.DataFrame(cn_amp_merged, columns=colnames)

    amp_regions['Chrom'] = amp_regions['Chrom'].astype(str)
    amp_regions['Start'] = pd.to_numeric(
        amp_regions['Start'], errors='coerce')
    amp_regions['End'] = pd.to_numeric(
        amp_regions['End'], errors='coerce')
    amp_regions = amp_regions.dropna(subset=['Start', 'End'])
    amp_regions['Start'] = amp_regions['Start'].astype(np.int64)
    amp_regions['End'] = amp_regions['End'].astype(np.int64)
    amp_regions.sort_values(colnames, kind='mergesort', inplace=True)
    amp_regions.reset_index(drop=True, inplace=True)

    amp_index = {}
    for chrom, region_df in amp_regions.groupby('Chrom', sort=False):
        amp_index[str(chrom)] = (
            region_df['Start'].to_numpy(dtype=np.int64, copy=False),
            region_df['End'].to_numpy(dtype=np.int64, copy=False),
        )

    amp_segment_set = set(
        amp_regions[colnames].itertuples(index=False, name=None))

    def component_ids(chrom, coords):
        """Return the containing or nearest boundary-adjacent CN component."""
        coords = np.asarray(coords, dtype=np.int64)
        out = np.full(coords.shape, -1, dtype=np.int64)
        region_data = amp_index.get(str(chrom))
        if region_data is None or coords.size == 0:
            return out

        starts, ends = region_data
        pos = np.searchsorted(starts, coords, side='right') - 1
        candidate_rows = np.flatnonzero(pos >= 0)
        if candidate_rows.size:
            candidate_pos = pos[candidate_rows]
            contained = ends[candidate_pos] >= coords[candidate_rows]
            out[candidate_rows[contained]] = candidate_pos[contained]

        outside_rows = np.flatnonzero(out < 0)
        for row in outside_rows:
            coord = coords[row]
            insert = np.searchsorted(starts, coord, side='right')
            best_component = -1
            best_distance = component_boundary_tolerance + 1

            previous = insert - 1
            if previous >= 0:
                distance = coord - ends[previous]
                if 0 < distance <= component_boundary_tolerance:
                    best_component = previous
                    best_distance = distance

            following = insert
            if following < starts.size:
                distance = starts[following] - coord
                if 0 < distance <= component_boundary_tolerance \
                        and distance < best_distance:
                    best_component = following

            out[row] = best_component
        return out

    # Keep only fields required for segment inference; avoid copying large
    # ReadID/Readbarcode payloads from breakpoint-pair tables.
    required_bp = ['Chrom', 'Coord', 'Clip', 'ClipDepth']
    missing = [x for x in required_bp if x not in bp_fine.columns]
    if missing:
        raise ValueError('bp_fine is missing columns: ' + ', '.join(missing))
    if restrict and 'CleanBP' not in bp_fine.columns:
        raise ValueError('bp_fine is missing CleanBP while restrict=True')

    bp_mask = bp_fine['CleanBP'].astype(bool) if restrict else slice(None)
    df = bp_fine.loc[bp_mask, required_bp].copy()
    df['Chrom'] = df['Chrom'].astype(str)
    df['Coord'] = pd.to_numeric(df['Coord'], errors='coerce')
    df['ClipDepth'] = pd.to_numeric(df['ClipDepth'], errors='coerce')
    df = df.dropna(subset=['Coord', 'ClipDepth'])
    df['Coord'] = df['Coord'].astype(np.int64)
    df['Clip'] = df['Clip'].astype(str)
    df.sort_values(['Chrom', 'Coord'], kind='mergesort', inplace=True)
    df.reset_index(drop=True, inplace=True)

    clip_breakpoint_set = set(
        df[['Chrom', 'Coord', 'Clip']].itertuples(
            index=False, name=None))

    pair_required = [
        'Chrom1', 'Coord1', 'Clip1', 'Chrom2', 'Coord2', 'Clip2']
    missing = [x for x in pair_required if x not in bp_pair.columns]
    if missing:
        raise ValueError('bp_pair is missing columns: ' + ', '.join(missing))
    if subseg and 'Count' not in bp_pair.columns:
        raise ValueError('bp_pair is missing Count while subseg=True')

    pair_columns = pair_required + [
        x for x in ('Count', 'Seq') if x in bp_pair.columns]
    pair_core = bp_pair.loc[:, pair_columns]

    if pair_core.empty:
        clip_pair = pair_core
    else:
        def pair_endpoint_is_clip(values):
            c1, p1, s1, c2, p2, s2 = values[:6]
            return (
                (str(c1), int(p1), str(s1)) in clip_breakpoint_set and
                (str(c2), int(p2), str(s2)) in clip_breakpoint_set)

        clip_mask = np.fromiter(
            (pair_endpoint_is_clip(row)
             for row in pair_core.itertuples(index=False, name=None)),
            dtype=bool, count=len(pair_core))
        clip_pair = pair_core.loc[clip_mask]

    if subseg:
        rm_bp_fine, inferred_records = get_segment_on_cnamp(
            amp_regions, clip_pair, knn=knn)
        rm_bp_fine_set = set(rm_bp_fine)
        if rm_amp_df and rm_bp_fine_set:
            keep = np.fromiter(
                (row not in rm_bp_fine_set
                 for row in df[['Chrom', 'Coord', 'Clip']].itertuples(
                     index=False, name=None)),
                dtype=bool, count=len(df))
            df = df.loc[keep].reset_index(drop=True)
    else:
        rm_bp_fine_set = set()
        inferred_records = []

    # Pair L/R boundaries only inside the same continuous CN component.
    for chrom, chrom_df in df.groupby('Chrom', sort=False):
        coords = chrom_df['Coord'].to_numpy(dtype=np.int64, copy=False)
        clips = chrom_df['Clip'].to_numpy(dtype=str, copy=False)
        depths = chrom_df['ClipDepth'].to_numpy(dtype=float, copy=False)
        components = component_ids(chrom, coords)

        for component in np.unique(components[components >= 0]):
            in_component = components == component
            component_coords = coords[in_component]
            component_clips = clips[in_component]
            component_depths = depths[in_component]

            left_mask = component_clips == 'L'
            right_mask = component_clips == 'R'
            if not left_mask.any() or not right_mask.any():
                continue

            left_coords = component_coords[left_mask]
            left_depths = component_depths[left_mask]
            right_coords = component_coords[right_mask]
            right_depths = component_depths[right_mask]

            for start in left_coords:
                key = (str(chrom), int(start), 'L')
                if key in rm_bp_fine_set:
                    continue
                first = np.searchsorted(
                    right_coords, start, side='right')
                stop = min(first + knn, right_coords.size)
                if first >= stop:
                    continue
                selected = int(np.argmax(right_depths[first:stop]))
                end = int(right_coords[first + selected])
                inferred_records.append([str(chrom), int(start), end])

            for end in right_coords:
                key = (str(chrom), int(end), 'R')
                if key in rm_bp_fine_set:
                    continue
                stop = np.searchsorted(left_coords, end, side='left')
                first = max(0, stop - knn)
                if first >= stop:
                    continue
                candidate_coords = left_coords[first:stop][::-1]
                candidate_depths = left_depths[first:stop][::-1]
                selected = int(np.argmax(candidate_depths))
                start = int(candidate_coords[selected])
                inferred_records.append([str(chrom), start, int(end)])

    # Direct loops already have two clip-derived endpoints, so scan only the
    # compact clip_pair table instead of the complete wide bp_pair table.
    n_pe_direct_loop_filtered = 0
    for row in clip_pair.itertuples(index=False):
        chrom1 = str(row.Chrom1)
        coord1 = int(row.Coord1)
        clip1 = str(row.Clip1)
        chrom2 = str(row.Chrom2)
        coord2 = int(row.Coord2)
        clip2 = str(row.Clip2)

        direct_seg = None
        if chrom1 == chrom2:
            if coord1 < coord2 and clip1 == 'L' and clip2 == 'R':
                direct_seg = [chrom1, coord1, coord2]
            elif coord1 > coord2 and clip1 == 'R' and clip2 == 'L':
                direct_seg = [chrom1, coord2, coord1]

        if direct_seg is None:
            continue
        is_pe_only = getattr(row, 'Seq', None) == 'PE_Support'
        direct_key = tuple(direct_seg)
        if is_pe_only and direct_key not in amp_segment_set:
            n_pe_direct_loop_filtered += 1
            continue
        inferred_records.append(direct_seg)

    if n_pe_direct_loop_filtered:
        print('filter PE-only direct loops without matching amplified '
              f'boundaries: {n_pe_direct_loop_filtered}')

    new_seg = pd.DataFrame.from_records(
        inferred_records, columns=colnames)
    new_seg['_type'] = 'inferred'

    if not new_seg.empty:
        new_seg['Chrom'] = new_seg['Chrom'].astype(str)
        new_seg['Start'] = pd.to_numeric(
            new_seg['Start'], errors='coerce')
        new_seg['End'] = pd.to_numeric(
            new_seg['End'], errors='coerce')
        new_seg = new_seg.dropna(subset=['Start', 'End'])
        new_seg['Start'] = new_seg['Start'].astype(np.int64)
        new_seg['End'] = new_seg['End'].astype(np.int64)

        segment_size = new_seg['End'] - new_seg['Start']
        new_seg = new_seg.loc[
            (segment_size >= min_segment_size) &
            (segment_size <= max_segment_size)].copy()

    # Defensive gate for subsegment/direct-loop records too. Measure the total
    # fraction of each inferred segment covered by merged amplified regions.
    if not new_seg.empty:
        component_supported = pd.Series(
            False, index=new_seg.index, dtype=bool)
        for chrom, row_index in new_seg.groupby(
                'Chrom', sort=False).groups.items():
            region_data = amp_index.get(str(chrom))
            if region_data is None:
                continue
            starts, ends = region_data
            rows = new_seg.loc[row_index]
            segment_starts = rows['Start'].to_numpy(
                dtype=np.int64, copy=False)
            segment_ends = rows['End'].to_numpy(
                dtype=np.int64, copy=False)
            segment_lengths = segment_ends - segment_starts
            amp_fraction = np.zeros(len(rows), dtype=float)

            for i, (segment_start, segment_end) in enumerate(zip(
                    segment_starts, segment_ends)):
                first = np.searchsorted(
                    ends, segment_start, side='right')
                stop = np.searchsorted(
                    starts, segment_end, side='left')
                if first >= stop or segment_lengths[i] <= 0:
                    continue
                overlaps = (
                    np.minimum(ends[first:stop], segment_end) -
                    np.maximum(starts[first:stop], segment_start))
                amp_fraction[i] = np.maximum(overlaps, 0).sum() / float(
                    segment_lengths[i])

            component_supported.loc[row_index] = (
                amp_fraction >= min_component_fraction)

        n_low_component_overlap = int((~component_supported).sum())
        if n_low_component_overlap:
            print('filter inferred segments below amplified-region overlap '
                  f'{min_component_fraction:.3f}: '
                  f'{n_low_component_overlap}')
        new_seg = new_seg.loc[component_supported].copy()

    old_seg = amp_regions.loc[:, colnames].copy()
    old_seg['_type'] = 'raw'

    segments = pd.concat(
        [old_seg, new_seg], ignore_index=True, copy=False)
    segments.drop_duplicates(
        subset=colnames, keep='first', inplace=True)
    segments.sort_values(
        colnames, kind='mergesort', inplace=True)
    segments.reset_index(drop=True, inplace=True)
    return segments

#------------------------------------------------------------------------------#
# get breakpoint from bp pairs and cn segments
#------------------------------------------------------------------------------#
def add_cn(df, cn_amp, bin_norm):
    """
    map copy number to breakpoint segments
    """

    # create GRange objects
    gr1 = GRange(df[['Chrom', 'Start', 'End']], 'dataframe_hasend')
    gr2 = GRange(cn_amp, 'dataframe_hasend')

    # # only use the mid point for the segments
    # op = []
    # for i in gr1.gr:
    #     mid = int((i[1].start + i[1].stop)/2)
    #     op.append((i[0], range(mid, mid+1), i[2]))
    # gr1.gr = op

    # extend binsize on both end of the cn_seg
    map_list = gr1.gmap(gr2, a_extend = 0, b_extend = 0, multi_hit = True)

    # output
    df['CN'] = [round(np.max([j[0] for j in i]))
        if (i is not None and i != []) else 0 for i in map_list]
    
    for row in df[df['CN'] == 0].itertuples():
        # get the bin_norm in seg
        tmp = bin_norm[(bin_norm['Chrom'] == row[1]) & 
                    (bin_norm['Coord'] >= row[2] - 10000) &
                    (bin_norm['Coord'] <= row[3] + 10000)]
        
        # mean
        cn_mean = np.round(tmp['CN'].mean()).astype('int')
        
        # fill the 0
        df.loc[row[0], 'CN'] = cn_mean

    return df


#------------------------------------------------------------------------------#
# annotate "bridge" segments whose interior is not amplified along its length
#------------------------------------------------------------------------------#
def filter_seg_by_interior_cn(
        seg, bin_norm,
        min_cn=5,
        min_amp_frac=0.5,
        return_all=False):
    """
    Annotate whether each segment is amplified along its interior.

    Adds:
        InteriorCNFrac:
            fraction of interior bins with CN >= min_cn
        InteriorCNPass:
            InteriorCNFrac >= min_amp_frac and segment CN > min_cn

    Only inferred segments (``_type == 'inferred'``) are evaluated. Raw
    segments are retained with InteriorCNFrac=1.0 and InteriorCNPass=True.
    Segment interval follows [Start, End), consistent with the original code.
    Inferred segments containing no bin are marked as passing.

    Set ``return_all=True`` to return every input row, including the
    ``InteriorCNFrac`` and ``InteriorCNPass`` annotations.  The default keeps
    the historical filtering behavior and output schema.
    """
    seg = seg.copy()

    if seg.empty:
        seg['InteriorCNFrac'] = pd.Series(index=seg.index, dtype=float)
        seg['InteriorCNPass'] = pd.Series(index=seg.index, dtype=bool)
        return seg

    if bin_norm.empty:
        seg['InteriorCNFrac'] = 1.0
        seg['InteriorCNPass'] = True
        return seg.reset_index(drop=True)

    required_seg = {'Chrom', 'Start', 'End'}
    required_bin = {'Chrom', 'Coord', 'CN'}

    if not required_seg.issubset(seg.columns):
        raise ValueError('seg must contain Chrom, Start and End')

    if not required_bin.issubset(bin_norm.columns):
        raise ValueError('bin_norm must contain Chrom, Coord and CN')

    # Avoid chromosome dtype mismatches.
    seg['Chrom'] = seg['Chrom'].astype(str)

    bins = bin_norm[['Chrom', 'Coord', 'CN']].copy()
    bins['Chrom'] = bins['Chrom'].astype(str)
    bins['Coord'] = pd.to_numeric(bins['Coord'], errors='coerce')
    bins['CN'] = pd.to_numeric(bins['CN'], errors='coerce')
    bins = bins.dropna(subset=['Coord'])
    bins['Coord'] = bins['Coord'].astype(np.int64)

    # Raw segments are trusted inputs and bypass the interior-CN filter. For
    # backward compatibility, old callers without _type are evaluated in full.
    seg['InteriorCNFrac'] = 1.0
    seg['InteriorCNPass'] = True
    if '_type' in seg.columns:
        inferred_mask = seg['_type'].astype(str).eq('inferred')
    else:
        inferred_mask = pd.Series(True, index=seg.index, dtype=bool)

    # Build one coordinate/prefix-sum index per chromosome.
    bin_index = {}
    for chrom, df in bins.groupby('Chrom', sort=False):
        df = df.sort_values('Coord')
        positions = df['Coord'].to_numpy(dtype=np.int64)
        amplified = ( df['CN'].to_numpy(dtype=float) >= min_cn).astype(np.int64)

        # prefix[i] = number of amplified bins before position i
        prefix = np.concatenate([ np.array([0], dtype=np.int64), np.cumsum(amplified, dtype=np.int64)])
        bin_index[chrom] = positions, prefix

    # Process inferred segments from one chromosome simultaneously.
    inferred_seg = seg.loc[inferred_mask]
    for chrom, row_index in inferred_seg.groupby('Chrom', sort=False).groups.items():
        data = bin_index.get(chrom)
        if data is None:
            continue
        positions, prefix = data
        rows = seg.loc[row_index]
        starts = pd.to_numeric(rows['Start'], errors='coerce').to_numpy(dtype=np.int64)
        ends = pd.to_numeric(rows['End'], errors='coerce').to_numpy(dtype=np.int64)

        # Select bins satisfying Start <= Coord < End.
        left = np.searchsorted(positions, starts, side='left')
        right = np.searchsorted(positions, ends, side='left')
        total_bins = right - left
        amplified_bins = (prefix[right] - prefix[left])

        # Preserve original behavior: no interior bin means pass.
        fractions = np.ones(total_bins.shape[0], dtype=float)
        has_bins = total_bins > 0
        fractions[has_bins] = (amplified_bins[has_bins] /total_bins[has_bins])
        seg.loc[row_index, 'InteriorCNFrac'] = fractions
        seg.loc[row_index, 'InteriorCNPass'] = fractions >= min_amp_frac

    # Apply the segment-level CN threshold only to inferred segments. Raw
    # segments keep the default pass state assigned above.
    segment_cn_pass = pd.to_numeric(seg['CN'], errors='coerce').gt(min_cn)
    seg['InteriorCNPass'] = (seg['InteriorCNPass'].astype(bool) & segment_cn_pass)
    seg['InteriorCNPass'] = seg['InteriorCNPass'].astype(bool)
    if return_all:
        return seg.reset_index(drop=True)
    # Keep the public segment table compatible with downstream code: the
    # filtering-only helper columns and _type are not returned.
    output_columns = ['Chrom', 'Start', 'End', 'CN']
    return seg.loc[seg['InteriorCNPass'], output_columns].copy().reset_index(drop=True)
