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
def bp_pair_filter(bp_duo, short_distance = 200, min_junc_cnt = 5,
                   max_offset_overhang = 1000, min_match_len = 5,
                   max_offset_insertion = 20, max_srp = 0.7):
    """
    filter breakpoint pairs

    max_srp : max simple repeat proportion
    """
    print(f'inital bp pair: {bp_duo.shape}')
    bp_pair = bp_duo[bp_duo['Count'] >= min_junc_cnt]
    print(f'inital bp pair: {bp_pair.shape}')

    ind = (bp_pair['Chrom1'] == bp_pair['Chrom2']) \
        & (abs(bp_pair['Coord2'] - bp_pair['Coord1']) < short_distance)
    print(f'remove based on short_distance: {sum(ind)}')

    # also control for max insertion size (+) and max overhang size (-)
    ind = ind | (bp_pair['offset'] > max_offset_insertion) \
              | (bp_pair['offset'] < -max_offset_overhang)
    print(f'remove based on max insertion size: {sum(ind)}')

    # get the matching sequence which in upper case
    match_seq = [re.search('([A-Z]+)', i).group(0) for i in bp_pair['Seq']]
    
    # length of the matched sequence
    match_len = [len(s) for s in match_seq]
    # get the whether is a simple repeat or not
    srp = [simple_repeat_proportion(s) for s in match_seq]
    # remove simple repeat bp pairs which is artificial
    ind2 = (np.array(srp) > max_srp) | \
           (np.array(match_len) < min_match_len)
    
    print(f'remove based on simple repeat: {sum(ind2)}')
    # include the PE-Support pairs
    ind = (~ind & (bp_pair['Seq'] == 'PE_Support')) | \
          (~ind & ~ind2)
    print(f'after filter: {bp_pair.loc[ind].shape}')
    return bp_pair.loc[ind].reset_index(drop=True)

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
def bp_refine(bp_pair, bp_cand_stats, cn_amp):
    """
    get breakpoints from breakpoint pairs and amplified segments
    """

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
    cn_amp_merged = misc.merge_bed(cn_amp)
    cn_amp_clip = dict()
    for i in cn_amp.itertuples():
        cn_amp_clip[(i[1], i[2], 'L')] = i[5]
        cn_amp_clip[(i[1], i[3], 'R')] = i[6]
        
    ## make breakpoints dataframe
    bplist = [i[1:4] for i in bp_stats.itertuples()]
    op = []
    for row in cn_amp_merged:
        L = (row[0], row[1], 'L')
        R = (row[0], row[2], 'R')
        if (L not in bplist):
            op.append([row[0], row[1], 'L', True, cn_amp_clip[L], 0, 0])
        if (R not in bplist):
            op.append([row[0], row[2], 'R', True, cn_amp_clip[R], 0, 0])
    cn_seg_df = pd.DataFrame(op, columns=bp_stats.columns)
    # cn_seg_df

    # merge and output
    df = pd.concat([bp_stats, cn_seg_df])
    df = df[pd.notnull(df.CleanBP)]
    df = df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
    return df.reset_index(drop=True), cn_amp_merged

#------------------------------------------------------------------------------#
# infer segments from breakpoint pairs
#------------------------------------------------------------------------------#
def get_segment_on_cnamp(cn_amp_merged, bp_pair_m, knn = 3):
    '''
    infer segments from bp pair table on cn ampilified

    usage:
        get_segment_on_cnamp(cn_amp, bp_pair)

    input:
    cn_amp: dataframe of cn segments
    bp_pair_m: pandas dataframe of bp pair
    '''
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
                if (bp_[1] == 55904441) or (bp_[1] == 56049369):
                    print(f'  ******* bp is L: {bp_}  from  {cand.iloc[ind[:knn],:]}*********')
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
                if (bp_[1] == 55904441) or (bp_[1] == 56049369):
                    print(f'  ******* bp is R: {bp_}  from  {cand.iloc[ind[:knn],:]}*********')
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

    for i in cn_amp_merged:
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

## 最新修改的地方
def get_segment(bp_fine, bp_pair, cn_amp_merged, restrict = False, \
        min_segment_size = 100, max_segment_size = 50000000, subseg = True, 
        rm_amp_df = True, knn = 3): # ******change
    """
    infer segments from breakpoint table

    usage:
        get_segment(bp_fine, bp_pair)

    input:
    bp_fine: dataframe of finally refined breakpoints
    output: pandas dataframe of segments
    restrict: only use breakpoints with True cleanBP to infer segments if True
    remove segments with size < min_segment_size and > max_segment_size
    """

    df = bp_fine.sort_values(['Chrom', 'Coord'])

    if (restrict): # only use CleanBP True
        df = df[df['CleanBP']]

    # init
    colnames = ['Chrom', "Start", "End"]

    if subseg == True:
        rm_bp_fine, seg = get_segment_on_cnamp(cn_amp_merged, bp_pair, knn = knn)  # ******change
        if rm_amp_df:
            filter_set = set(rm_bp_fine)
            def filter_df(filter_set, df):
                df1 = df.set_index(['Chrom', 'Coord', 'Clip'])
                df = df[~df1.index.isin(filter_set)].reset_index(drop=True)
                return df
            df = filter_df(filter_set, df)
    else: 
        rm_bp_fine = ['Chrom']
        seg = []

    # infer potential segments from breakend patterns (<-...->)
    for chrom in set(df['Chrom']):
        ind = df['Chrom'] == chrom
        df_all = df.loc[ind].copy()
        # df_all.loc[df_all['InDepth'] == 0, 'CleanBP'] = False
        df_l = df_all.loc[df_all['Clip'] == "L"].reset_index(drop=True)
        df_r = df_all.loc[df_all['Clip'] == "R"].reset_index(drop=True)
    
        if (df_l.size > 0 and df_r.size > 0):
            for d,row in df_l.iterrows():
                if tuple(row[0:3]) not in rm_bp_fine:
                    dst = (row['Coord']-df_r['Coord'])
                    # dst = np.where(dst < 0, -dst, np.inf)
                    cand = df_r.loc[dst < 0].reset_index(drop=True)
                    sort_index = np.argsort(cand['Coord'])
                    # sort_index = np.argsort(dst)
                    try:
                        cand = cand.iloc[sort_index[:knn],:].reset_index(drop=True)
                        # ind = list(df_r['Count'][sort_index[:5]])
                        ind = np.argsort(-cand['ClipDepth'])[0]
                        # end = [cand['Coord'][ind] for i in range(ind+1)]
                        end = cand['Coord'][ind]
                        seg.append(list(row[0:2]) + [end])
                        # for i in end:
                        #     seg.append(list(row[0:2]) + [i])
                    except:
                        1
                        # print(row, '\n')

            for d,row in df_r.iterrows():
                if tuple(row[0:3]) not in rm_bp_fine:
                    dst = (df_l['Coord']-row['Coord'])
                    # dst = np.where(dst < 0, -dst, np.inf)
                    cand = df_l.loc[dst < 0].reset_index(drop=True)
                    sort_index = np.argsort(-cand['Coord'])
                    try:
                        cand = cand.iloc[sort_index[:knn], :].reset_index(drop=True)
                        ind = np.argsort(-cand['ClipDepth'])[0]
                        end = cand['Coord'][ind]
                        seg.append([row[0]] + [end] + [row[1]])
                        # ind = list(df_l['CleanBP'][sort_index]).index(True)
                        # end = [df_l['Coord'][sort_index[i]] for i in range(ind+1)]
                        # for i in end:
                        #    seg.append([row[0]] + [i] + [row[1]])
                    except:
                        1
                        # print(row, '\n')

    # add segment from direct loop
    # in case detailed structure can't be found
    # return this simple circle
    for row in bp_pair.itertuples():
        if (row[1] == row[4]):
            if (row[2] < row[5] and row[3] == 'L' and row[6] == 'R'):
                seg.append([row[1], row[2], row[5]])
            elif (row[2] > row[5] and row[3] == 'R' and row[6] == 'L'):
                seg.append([row[1], row[5], row[2]])

    # merge all segments
    seg = pd.DataFrame(seg, columns = colnames).drop_duplicates()

    # remove oversized segments
    seg_size = seg['End'] - seg['Start']
    seg = seg[(seg_size >= min_segment_size) & (seg_size <= max_segment_size)]

    return seg.sort_values(['Chrom', 'Start', 'End']).reset_index(drop=True)

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
