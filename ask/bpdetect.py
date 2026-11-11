################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd
from itertools import combinations
#------------------------------------------------------------------------------#
from grange import GRange
from covbam import check_read

#------------------------------------------------------------------------------#
# filter breakpoint by "n" supporting reads
#------------------------------------------------------------------------------#
def bp_filter_by_n(bpall, n = 5, f_sort = True, pe_window = 150, edges = None,  min_edge_support = 2):
    """
    filter breakpoint by # of softclip and hardclip reads
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    usage:
        breakpoint_filter_by_n(bpall)

    input: output of breakpoint_from_bam
    output: pandas dataframe of split read counts in a dict

    CHANGE (20260718, ported from ask3-EGA 20260705):
        Keep breakpoints supported by at least n clip reads. When precomputed
        discordant edges are provided, also retain a breakpoint with 1 to n-1
        clip reads if a compatible edge endpoint is found on the same
        chromosome, with the same orientation, within pe_window, and the edge
        has at least min_edge_support supporting reads.

        This preserves the base-pair-level coordinate supplied by clip reads
        while allowing discordant PE evidence to rescue low-clip junctions.
        Endpoint coordinates are indexed in sorted arrays and matched with
        searchsorted, avoiding repeated BAM scans. If edges is None, the
        original clip-count threshold (Count >= n) is used.
    """
    df_left = bpall['left']
    df_right = bpall['right']

    if (n > 0):
        if edges is not None:
            # In an edge, '+' denotes a right endpoint (R) and '-' denotes a
            # left endpoint (L). Build the endpoint index once.
            edge_index = {}
            for edge in edges:
                if len(edge) < 3 or edge[2] < min_edge_support:
                    continue
                for chrom, coord, strand in edge[:2]:
                    clip = 'R' if strand == '+' else 'L'
                    edge_index.setdefault((chrom, clip), []).append(int(coord))
            edge_index = {key: np.sort(np.asarray(pos, dtype=np.int64))
                          for key, pos in edge_index.items()}

            def _keep_by_nearby_edge(df):
                """Keep Count >= n or rescue a low-count clip near a compatible edge."""
                if df.empty:
                    return df
                out = df.copy().reset_index(drop=True)
                keep = out['Count'].ge(n)
                low = out.loc[out['Count'].ge(1) & out['Count'].lt(n)]
                rescued = 0
                for (chrom, clip), group in low.groupby(['Chrom', 'Clip'], sort=False):
                    pos = edge_index.get((chrom, clip))
                    if pos is None or pos.size == 0:
                        continue
                    coord = group['Coord'].to_numpy(dtype=np.int64)
                    ins = np.searchsorted(pos, coord)
                    dist = np.full(coord.size, np.iinfo(np.int64).max, dtype=np.int64)
                    valid = ins < pos.size
                    dist[valid] = np.abs(pos[ins[valid]] - coord[valid])
                    valid = ins > 0
                    dist[valid] = np.minimum(
                        dist[valid], np.abs(pos[ins[valid] - 1] - coord[valid]))
                    rescued_idx = group.index[dist <= pe_window]
                    keep.loc[rescued_idx] = True
                    rescued += rescued_idx.size
                print(f'  [bp_filter_by_n] clip>=n retained '
                      f'{int(out["Count"].ge(n).sum())}; edge({pe_window} bp) '
                      f'rescued {rescued}/{low.shape[0]}')
                return out.loc[keep]

            df_left = _keep_by_nearby_edge(df_left)
            df_right = _keep_by_nearby_edge(df_right)
        else:
            df_left = df_left[df_left['Count'] >= n]
            df_right = df_right[df_right['Count'] >= n]

    if (f_sort):
        df_left = df_left.sort_values('Count', ascending=False)
        df_right = df_right.sort_values('Count', ascending=False)

    ## return a dict of chrom.position with "left" and "right" keys
    return dict(left=df_left, right=df_right)


# #------------------------------------------------------------------------------#
# # CHANGE(20260718-NEW): Helpers for joint clip and discordant-PE filtering
# # ported from ask3-EGA
# #------------------------------------------------------------------------------#
# def _build_amp_lookup(cn_amp):
#     """Convert cn_amp to {chrom: [(start, end), ...]} for interval lookup."""
#     amp = {}
#     if cn_amp is None or len(cn_amp) == 0:
#         return None
#     for row in cn_amp.itertuples():
#         chrom, start, end = row[1], int(row[2]), int(row[3])
#         amp.setdefault(chrom, []).append((start, end))
#     return amp


# def _in_amp(amp, chrom, pos, flank):
#     """Check whether pos +/- flank overlaps an amplified segment."""
#     if amp is None:
#         return True
#     for s, e in amp.get(chrom, ()):
#         if pos >= s - flank and pos <= e + flank:
#             return True
#     return False


# def _build_boundary_lookup(cn_amp):
#     """Build {chrom: sorted[boundary_pos, ...]} for CN segment boundaries.
#     Rank PE rescue candidates by distance to the nearest boundary."""
#     bnd = {}
#     if cn_amp is None or len(cn_amp) == 0:
#         return None
#     for row in cn_amp.itertuples():
#         chrom, start, end = row[1], int(row[2]), int(row[3])
#         bnd.setdefault(chrom, []).extend((start, end))
#     for c in bnd:
#         bnd[c] = sorted(bnd[c])
#     return bnd


# def _dist_to_boundary(bnd, chrom, pos):
#     """Return the distance from pos to the nearest CN boundary."""
#     import bisect
#     if bnd is None:
#         return 0
#     arr = bnd.get(chrom)
#     if not arr:
#         return 1 << 30
#     i = bisect.bisect_left(arr, pos)
#     best = 1 << 30
#     if i < len(arr):
#         best = min(best, abs(arr[i] - pos))
#     if i > 0:
#         best = min(best, abs(pos - arr[i - 1]))
#     return best


# def _count_discordant_pe_at(bamf, chrom, pos, clip, window, min_mate_dist):
#     """
#     Count discordant PE reads supporting breakpoint (chrom, pos, clip).
#     Follow the same convention as bpjoint.discover_bp_from_discordant:
#       clip=='R': a forward read with reference_end near pos;
#       clip=='L': a reverse read with reference_start near pos.
#     A pair is discordant when it is not proper and its mate is on another
#     chromosome or farther than min_mate_dist on the same chromosome.
#     """
#     cnt = 0
#     for read in bamf.fetch(chrom, max(0, pos - window), pos + window):
#         if not check_read(read):
#             continue
#         if read.is_proper_pair or read.mate_is_unmapped \
#                 or read.next_reference_name is None:
#             continue
#         far = (read.next_reference_name != chrom) or \
#               (abs(read.next_reference_start - read.reference_start) > min_mate_dist)
#         if not far:
#             continue
#         if clip == 'R':
#             if (not read.is_reverse) and abs(read.reference_end - pos) <= window:
#                 cnt += 1
#         else:  # 'L'
#             if read.is_reverse and abs(read.reference_start - pos) <= window:
#                 cnt += 1
#     return cnt


# def _keep_by_clip_plus_pe(df, n, bamfile, cn_amp, pe_window, min_mate_dist, max_pe_cand):
#     """
#     Keep breakpoints with clip>=1 and clip + discordant_PE >= n.
#     Keep Count>=n directly. Rescue 1<=Count<n candidates with PE evidence
#     only inside amplified regions and within max_pe_cand.
#     """
#     if df.empty:
#         return df
#     pass_clip = df[df['Count'] >= n]          # clip support passes directly
#     need_pe = df[df['Count'] < n]             # requires PE support
#     if need_pe.empty:
#         return pass_clip

#     amp = _build_amp_lookup(cn_amp)
#     # Count PE support only inside amplified regions.
#     if amp is not None:
#         mask = [_in_amp(amp, r.Chrom, int(r.Coord), pe_window)
#                 for r in need_pe.itertuples()]
#         need_pe = need_pe[mask]
#     if need_pe.empty:
#         return pass_clip

#     # Limit large candidate sets. Rank by distance to the nearest CN boundary,
#     # then by clip support, because true junctions occur at segment boundaries.
#     if need_pe.shape[0] > max_pe_cand:
#         print(f'  [bp_filter_by_n] {need_pe.shape[0]} PE rescue candidates; '
#               f'keep the closest {max_pe_cand} to CN boundaries')
#         bnd = _build_boundary_lookup(cn_amp)
#         need_pe = need_pe.assign(
#             _bdist=[_dist_to_boundary(bnd, r.Chrom, int(r.Coord))
#                     for r in need_pe.itertuples()]
#         ).sort_values(['_bdist', 'Count'], ascending=[True, False]).head(max_pe_cand)
#         need_pe = need_pe.drop(columns='_bdist')

#     rescued_idx = []
#     with pysam.AlignmentFile(bamfile, 'rb') as bamf:
#         for r in need_pe.itertuples():
#             pe = _count_discordant_pe_at(
#                 bamf, r.Chrom, int(r.Coord), r.Clip, pe_window, min_mate_dist)
#             if r.Count + pe >= n:
#                 rescued_idx.append(r.Index)
#     rescued = need_pe.loc[rescued_idx]
#     print(f'  [bp_filter_by_n] clip>=n kept {pass_clip.shape[0]}; '
#           f'PE rescued {rescued.shape[0]}/{need_pe.shape[0]}')
#     return pd.concat([pass_clip, rescued]) if not rescued.empty else pass_clip


#------------------------------------------------------------------------------#
# at most include 200 bp candidates on each side
#------------------------------------------------------------------------------#

def top_n_scores(group,byname,top_n):
    '''
    take top N by group
    '''
    return group.nlargest(top_n, byname)


def stratified_sampling(data, group_column, value_column, sample_percentage=0.3,even=False,seed=1):
    """
    stratified sampling by group

    Parameters:
    data (DataFrame): dataframe
    group_column (str): sample by column
    value_column (str): to be computed proportion 
    sample_percentage (float): sample percentage
    even: wheater the same number for each group
    seed (int, optional): 

    Return:
    sample (DataFrame): 
    """
    np.random.seed(seed)
    
    # caculated the proportion of each group
    probabilities = data[value_column] / data.groupby(group_column)[value_column].transform('sum')
    probabilities = probabilities / probabilities.sum()
    
    # the sample size of each group
    # if even:
    #     chr_num = Counter(data[group_column])
    #     percentage_by_chr = sample_percentage / len(chr_num)
    #     group_sizes = [int(percentage_by_chr * count) for count in chr_num.values()]
    # else:
    #     chr_num = Counter(data[group_column])
    #     tmp = sample_percentage * data.shape[0] / len(chr_num)  
    #     group_sizes = [ int(tmp) + 1 for _ in range(len(chr_num))]
    #     del tmp
    
    group_sizes = int(sample_percentage * data.shape[0])


    # sample with numpy.random.choice
    sampled_indices = np.random.choice(data.index, size=group_sizes, replace=False, p=probabilities)

    sampled_data = data.loc[sampled_indices]
    
    return sampled_data


def stratified_sampling_by_intervals(data, value_column, num_intervals, 
                                      num_samples_per_interval,seed=1):
    """
    Parameter:
    data (DataFrame):data
    value_column (str): the column used to split interval.
    num_samples_per_interval (int): the number of samples for each interval
    seed (int, optional): random seed

    Return:
    sample (DataFrame): data after sample.
    """
    np.random.seed(seed)
    
    # sort by value_column
    sorted_df = data.sort_values(by=value_column)
    
    # compute the length of interval
    interval_length = len(sorted_df) // num_intervals
    
    # sample in each interval
    samples = []
    for i in range(num_intervals):
        start = i * interval_length
        end = (i + 1) * interval_length
        interval_data = sorted_df.iloc[start:end]
        if num_samples_per_interval >= interval_data.shape[0]:
            num_samples_per_interval = interval_data.shape[0]
        interval_sample = interval_data.sample(n=num_samples_per_interval, random_state=seed, replace=False)
        samples.append(interval_sample)
    
    sub_data = pd.concat(samples).reset_index(drop=True)
    return sub_data

def bp_top_n(bp_cand, topn = 200, num_chrom = 200, num_stratify = None):
    """
    at most include 200 bp candidates on each side
    """
    print(f'bp_top_n: 1 round')
    res = []
    left = bp_cand['left']
    left['Group'] = 1
    if (len(left.index) > topn):
        cutoff = left['Count'][topn-1]
        res.append(left[left['Count'] > cutoff])
    else:
        res.append(left)

    right = bp_cand['right']
    right['Group'] = 1
    if (len(right.index) > topn):
        cutoff = right['Count'][topn-1]
        res.append(right[right['Count'] > cutoff])
    else:
        res.append(right)

    # sample by chrom
    if num_chrom:
        if (len(left.index) > topn) or (len(right.index) > topn):
            print(f'bp_top_n: 2 round')
            if (len(left.index) > topn):
                df1 = bp_cand['left'].groupby('Chrom', \
                                        group_keys=False, sort=False).apply(top_n_scores, 'Count', topn)
            else:
                df1 = pd.DataFrame()

            if (len(right.index) > topn):
                df2 = bp_cand['right'].groupby('Chrom', \
                                        group_keys=False, sort=False).apply(top_n_scores, 'Count', topn)
            else:
                df2 = pd.DataFrame()

            df1['Group'] = 2
            df2['Group'] = 2
            res.append(df1)
            res.append(df2)

    if num_stratify:
        if (len(left.index) > topn) or (len(right.index) > topn):
            print(f'bp_top_n: 3 round')
            # by intervals
            df1 = stratified_sampling_by_intervals(left,'Count', 3, num_stratify)
            df2 = stratified_sampling_by_intervals(right,'Count', 3, num_stratify)
            df1['Group'] = 3
            df2['Group'] = 3
            res.append(df1)
            res.append(df2)

    columns =['Chrom','Coord','Clip','Count','Group']
    df = pd.concat(res).sort_values(columns[:3]).reset_index(drop=True)
    return df     # dict(left = left,right = right)

# def bp_top_n(bp_cand, topn = 200):
#     """
#     at most include 200 bp candidates on each side
#     """
#     left = bp_cand['left']
#     if (len(left.index) > topn):
#         cutoff = left['Count'][topn-1]
#         left = left[left['Count'] > cutoff]

#     right = bp_cand['right']
#     if (len(right.index) > topn):
#         cutoff = right['Count'][topn-1]
#         right = right[right['Count'] > cutoff]

#     return dict(left=left, right=right)

#------------------------------------------------------------------------------#
# filter breakpoints by blacklist intervals
#------------------------------------------------------------------------------#
def bp_filter_by_blacklist(df, blacklistfile, f_sort=True, isdf = True):
    """
    filter breakpoint by blacklist intervals
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    usage:
        breakpoint_filter_by_blacklist(bpall)

    input: output of breakpoint_from_bam
    output: pandas dataframe of split read counts in a dict

    """
    if isdf:
        gr1 = GRange(df, 'bpdataframe')
        gr2 = GRange(blacklistfile, 'bedfile')
        # extend binsize bp on both end of the cn_amplicon
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 100, invert = True, mode=2)

        # a = pd.DataFrame(\
        #     [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
        #     , columns = df.columns)
        df = df[gr]
        if (f_sort):
            df = df.sort_values('Count', ascending=False)
        return df
    else:
        def filter_blacklist(df, blacklistfile):
            gr1 = GRange(df, 'dataframe')
            gr2 = GRange(blacklistfile, 'bedfile')
            # extend 100 bp on both end of the blacklist
            gr = gr1.intersect(gr2, a_extend = 0, b_extend = 100, invert = True, mode=1)
            return pd.DataFrame(\
                [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
                , columns=df.columns)

        df_left = df['left']   # df is bp_all
        df_right = df['right']

        df_left = filter_blacklist(df_left, blacklistfile)
        df_right = filter_blacklist(df_right, blacklistfile)

        if (f_sort):
            df_left = df_left.sort_values('Count', ascending=False)
            df_right = df_right.sort_values('Count', ascending=False)

        ## return a dict of chrom.position with "left" and "right" keys
        return dict(left=df_left, right=df_right)


#------------------------------------------------------------------------------#
# filter breakpoints by amplified segments
#------------------------------------------------------------------------------#
def bp_filter_by_amplicon(df, cn_amplicon, binsize = 10000, f_sort = True, isdf = True):  #******change*****
    """
    filter breakpoint by amplified segments
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    bin_extend: extend cn_amplicon by <N> bp on both sides
    """

    if isdf:
        gr1 = GRange(df, 'bpdataframe')
        gr2 = GRange(cn_amplicon, 'dataframe_hasend')
        # extend binsize bp on both end of the cn_amplicon
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 2*binsize, invert = False, mode=2)

        # a = pd.DataFrame(\
        #     [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
        #     , columns = df.columns)
        df = df[gr]

        if (f_sort):
            df = df.sort_values('Count', ascending=False)
        return df
    else:
        def filter_amplicon(df, cn_amplicon):
            gr1 = GRange(df, 'dataframe')
            gr2 = GRange(cn_amplicon, 'dataframe_hasend')
            # extend binsize bp on both end of the cn_amplicon
            gr = gr1.intersect(gr2, a_extend = 0, b_extend = 2*binsize, invert = False, mode = 1)
            gr.gr
            return pd.DataFrame(\
                [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
                , columns=df.columns)
        df_left = df['left']   # bp_all
        df_right = df['right']

        df_left = filter_amplicon(df_left, cn_amplicon)
        df_right = filter_amplicon(df_right, cn_amplicon)

        if (f_sort):
            df_left = df_left.sort_values('Count', ascending=False)
            df_right = df_right.sort_values('Count', ascending=False)

        ## return a dict of chrom.position with "left" and "right" keys
        return dict(left=df_left, right=df_right)


#------------------------------------------------------------------------------#
# Breakpoint evaluation by sequence depth
#------------------------------------------------------------------------------#
def bp_seq_depth(bp_cand_df, bamfile, perc = 0.2, mapq = 20, nmmax = 1, isdf =  True):
    """
    calculate the sequence depth on the breakpoint (InDepth)
    and 1-bp out of breakpoint (OutDepth)

    Clean breakpoint (CleanBP) is defined as the
    ClipDepth/InDepth >= perc,
    and means it can't be set as alternative segment boundary

    """    
    # bp_duo is bp_cand from bp_all
    bp_cand_df_ = bp_cand_df.drop_duplicates(['Chrom', 'Coord', 'Clip'])
    
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        op = []
        for index, row in bp_cand_df_.iterrows():  # here is can change
            # depth_in = bamf.count(row[0], row[1], row[1]+1, read_callback = check_read)

            depth_in = len([
                read for read in bamf.fetch(row[0], row[1], row[1]+1) \
                if check_read(read, mapq = mapq, nmmax = nmmax)])

            depth_out = depth_in - row[3]
            # if (depth_out < 0):
            #     print([[row], depth_in])

            depth_bp = [depth_in, depth_out, row[3]/depth_in >= perc]
            
            op.append(depth_bp)
    colnames=['InDepth', 'OutDepth', 'CleanBP']
    df = pd.concat([bp_cand_df_.reset_index(drop=True), pd.DataFrame(op, columns=colnames)], axis=1)
    df = pd.merge(bp_cand_df, df, on = ['Chrom', 'Coord', 'Clip', 'Count'],  how = 'left').reset_index(drop=True)
    df = df[['Chrom', 'Coord', 'Clip', 'CleanBP', 'Count', 'InDepth', 'OutDepth', 'Group_x']]
    df.rename(columns={'Count':'ClipDepth', 'Group_x': 'Group'}, inplace=True)
    return df

#------------------------------------------------------------------------------#
# filter bp cand based on  bp from bam
#------------------------------------------------------------------------------#
# def filter_bp_cand(bp_cand_df, bp_cand_bam):
#     bp_cand_df = bp_cand_df.sort_values('Count', ascending=False).reset_index(drop=True)
#     bp_cand_df_ = pd.merge(bp_cand_df, 
#                             bp_cand_bam[['Chrom','Coord','Clip']], on = ['Chrom','Coord','Clip'],\
#                             indicator = True, how = 'left')\
#                                 .query('_merge == "left_only"').drop('_merge', axis = 1).reset_index(drop=True)
#     return bp_cand_df_

def refine_bp_cand_bam(bp_duo, bp_all):
    '''
    refine the Count of bp clip by bam based on bp_all
    '''
    col1 = ['Chrom1', 'Coord1', 'Clip1', 'Count']
    col2 = ['Chrom2', 'Coord2', 'Clip2','Count']
    newcol = ['Chrom', 'Coord', 'Clip', 'Count']
    df = pd.concat([bp_duo[col1].set_axis(newcol, axis=1), \
                    bp_duo[col2].set_axis(newcol, axis=1)], copy=False)\
                        .groupby(newcol[:3], as_index = False).agg({'Count': ['sum']}).set_axis(newcol, axis=1)
    df  = df.set_index(['Chrom','Coord','Clip']) 

    bp_all_df = pd.concat([bp_all['left'][bp_all['left']['Count']>1], \
                           bp_all['right'][bp_all['right']['Count']>1]], copy=False).set_index(['Chrom','Coord','Clip'])

    # merge blocks
    df = pd.merge(df, bp_all_df, left_index=True, right_index=True, how='left', copy=False)

    # fill missing values
    df['Count_y'] = df['Count_y'].fillna(df['Count_x'])

    # take the maximum count
    df['Count'] = np.maximum(df['Count_x'], df['Count_y'])

    # drop intermediate columns
    df = df.drop(['Count_x', 'Count_y'], axis=1).reset_index()

    df['CleanBP'] =  True
    df['OutDepth'] =  0
    df['InDepth'] =  0
    df = df[['Chrom', 'Coord', 'Clip', 'CleanBP', 'Count', 'InDepth', 'OutDepth']]
    df.rename(columns = {'Count':'ClipDepth'}, inplace=True)
    return df
