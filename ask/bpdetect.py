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
def bp_filter_by_n(bpall, n=5, f_sort=True):
    """
    filter breakpoint by # of softclip and hardclip reads
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    usage:
        breakpoint_filter_by_n(bpall)

    input: output of breakpoint_from_bam
    output: pandas dataframe of split read counts in a dict

    """
    df_left = bpall['left']
    df_right = bpall['right']

    if (n > 0):
        df_left = df_left[df_left['Count']>=n]
        df_right = df_right[df_right['Count']>=n]
    if (f_sort):
        df_left = df_left.sort_values('Count', ascending=False)
        df_right = df_right.sort_values('Count', ascending=False)

    ## return a dict of chrom.position with "left" and "right" keys
    return dict(left=df_left, right=df_right)


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


    # 使用 numpy.random.choice 进行抽样
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

    # 分块合并
    df = pd.merge(df, bp_all_df, left_index=True, right_index=True, how='left', copy=False)

    # 使用map填充  
    df['Count_y'] = df['Count_y'].fillna(df['Count_x'])

    # 取最大count
    df['Count'] = np.maximum(df['Count_x'], df['Count_y'])

    # 删除中间列
    df = df.drop(['Count_x', 'Count_y'], axis=1).reset_index()

    df['CleanBP'] =  True
    df['OutDepth'] =  0
    df['InDepth'] =  0
    df = df[['Chrom', 'Coord', 'Clip', 'CleanBP', 'Count', 'InDepth', 'OutDepth']]
    df.rename(columns = {'Count':'ClipDepth'}, inplace=True)
    return df