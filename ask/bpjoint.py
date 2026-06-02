################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd
import os
import re
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from itertools import combinations
from multiprocessing.pool import ThreadPool as Pool
from functools import partial
#------------------------------------------------------------------------------#
from covbam import check_read
import misc
import logger

logg = logger.Logger()
#------------------------------------------------------------------------------#
def check_read_clip(read, clip_type = 'soft_left', seq = False):
    """
    true if the read is the clip read
    """
    if seq:
        tf = check_read(read, mapq=0)
    else:
        tf = check_read(read)

    if (tf):
        clip_tf = False
        if (clip_type == 'soft_left' and read.cigartuples is not None):
            clip_tf = read.cigartuples[0][0] == 4 or (read.cigartuples[-1][0] == 5 and seq) or (read.cigartuples[0][0] == 5 and seq)
        elif (clip_type == 'soft_right' and read.cigartuples is not None):
            clip_tf = read.cigartuples[-1][0] == 4 or (read.cigartuples[-1][0] == 5 and seq)
        tf = tf and clip_tf

    return tf

#------------------------------------------------------------------------------#
def get_consensus_sequence(bamfile, contig, start, stop):
    """
    get the consensus sequence of a given region
    """
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        seq = bamf.count_coverage(contig, start, stop)
        nn = ['A', 'C', 'G', 'T']
        return ''.join([nn[i] for i in np.argmax(seq, axis=0)])

#------------------------------------------------------------------------------#
def rev_compl(seq):
    nn = defaultdict(
        lambda: 'N', {
            'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
            'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'})
    return ''.join([nn[a] for a in seq][::-1])


#------------------------------------------------------------------------------#
# main function to joint any two breakpoints
#------------------------------------------------------------------------------#
def join_breakpoint(bamfile, a, b, rm_clipped_name = None, min_nt = 2, seq_len = 200, \
    offset = 0, match_method = 'fuzzy_match'):
    """
    detect whether the two breakpoints are joined in the genome

    a, b - breakpoints in a tuple, such as:
        a = ('chr7', 54830975, 'L')
        b = ('chr7', 56117062, 'R')

    pattern -  mean the clip patterns of a and b breakpoints
        'LR': a - left clip; b - right clip
        'LL': a - left clip; b - left clip
        'RR': a - right clip; b - right clip

    min_qlen - minimal length of clipped bases to compare
        remove read with <= min_qlen

    seq_len - # of bases within the breakpoints to compare

    match_method - perfect_match, partial_match, fuzzy_match

    usage:
        join_breakpoint(bamfile, ('chr7', 55194959, 'L'), ('chr7', 56117062, 'R'))
    """

    if (a[2] == 'L'):
        
        a_seq, a_clipped, a_clipped_name = left_clip(bamfile, a, seq_len)
    else:
        a_seq, a_clipped, a_clipped_name = right_clip(bamfile, a, seq_len)

    if (b[2] == 'L'):
        b_seq, b_clipped, b_clipped_name = left_clip(bamfile, b, seq_len)
    else:
        b_seq, b_clipped, b_clipped_name = right_clip(bamfile, b, seq_len)
    # with open('/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_1113/ask3_both_knn3_mr10_GBM39_1113_nofilter_circ/reads_55222714.txt', 'a') as f:
    #     f.write(f'{a}\t{b}\t{a_clipped_name}\t{b_clipped_name}\n')
    # remove too short reads
    if (min_nt > 0):
        if rm_clipped_name:
            def filter_clipped(clipped, clipped_name, min_nt):
                clipped_ = []
                clipped_name_ = []
                for i, read in enumerate(clipped):
                    if (len(read) >= min_nt):
                        clipped_.append(read)
                        clipped_name_.append(clipped_name[i])
                return clipped_, clipped_name_
            a_clipped, a_clipped_name = filter_clipped(a_clipped, a_clipped_name, min_nt)
            b_clipped, b_clipped_name = filter_clipped(b_clipped, b_clipped_name, min_nt)
        else:
            a_clipped = [i for i in a_clipped if len(i) >= min_nt]
            b_clipped = [i for i in b_clipped if len(i) >= min_nt]
            # a_clipped_name = [a_clipped_name[i] for i in range(len(a_clipped)) \
            #                   if len(a_clipped[i]) >= min_nt]
            # b_clipped_name = [b_clipped_name[i] for i in range(len(b_clipped)) \
            #             if len(b_clipped[i]) >= min_nt]
    
    if rm_clipped_name:
        def filter_clip_by_bpp(rm_clipped_name, clipped_name, clipped):
            clipped_ = []
            clipped_name_ = []
            num = 0
            for i, name in enumerate(clipped_name):
                if name in rm_clipped_name:
                    num += 1
                else:
                    clipped_.append(clipped[i])    
                    clipped_name_.append(name)
            return clipped_, clipped_name_, num
        #  remove the clip read that exists in bp pairs from bamfile
        a_clipped, a_clipped_name, num1 = filter_clip_by_bpp(rm_clipped_name, a_clipped_name, a_clipped)
        b_clipped, b_clipped_name, num2 = filter_clip_by_bpp(rm_clipped_name, b_clipped_name, b_clipped)
    else:
        num1 = 0
        num2 = 0

    # format clip reads and breakend sequence to the same format
    if (a[2]+b[2] == 'LR'):
        _b2a = R2L(b_clipped, a_seq) # map b clipped to  a-seq
        _a2b = L2R(a_clipped, b_seq)
    elif (a[2]+b[2] == 'RL'):
        _b2a = L2R(b_clipped, a_seq)
        _a2b = R2L(a_clipped, b_seq)
    elif (a[2]+b[2] == 'LL'):
        _b2a = L2L(b_clipped, a_seq)
        _a2b = L2L(a_clipped, b_seq)
    elif (a[2]+b[2] == 'RR'):
        _b2a = R2R(b_clipped, a_seq)
        _a2b = R2R(a_clipped, b_seq)
    else:
        raise(Exception("Invalid pattern"))

    # map clip reads to breakend sequence
    if (match_method == 'perfect_match'):
        b2a = perfect_match(_b2a[0], _b2a[1], num2, offset = offset, min_nt = min_nt)
        a2b = perfect_match(_a2b[0], _a2b[1], num1, offset = offset, min_nt = min_nt)
    elif (match_method == 'partial_match'):
        b2a = partial_match(_b2a[0], _b2a[1], num2, min_nt = min_nt)
        a2b = partial_match(_a2b[0], _a2b[1], num1, min_nt = min_nt)
    elif (match_method == 'fuzzy_match'):
        b2a = fuzzy_match(_b2a[0], _b2a[1], num2, offset = offset, min_nt = min_nt)
        a2b = fuzzy_match(_a2b[0], _a2b[1], num1, offset = offset, min_nt = min_nt)
    else:
        raise(Exception("Invalid method"))

    return [b2a, a2b]
#------------------------------------------------------------------------------#
def left_clip(bamfile, l, seq_len):   # ************change
    # get the clip read for breakpoint 
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        l_read = [
            read for read in bamf.fetch(l[0], l[1], l[1]+1)
            if read.reference_start == l[1]
            and check_read_clip(read, 'soft_left')]

    # get the sequence bases within the breakpoint
    l_seq = get_consensus_sequence(bamfile, l[0], l[1], l[1]+seq_len)

    # get the clipped reads of the breakpoint
    l_clipped = [read.query_sequence[0:read.query_alignment_start]
                 for read in l_read]
    
    l_clipped_name = [read.query_name for read in l_read]

    return [l_seq, l_clipped, l_clipped_name]

#------------------------------------------------------------------------------#
def right_clip(bamfile, r, seq_len):   # ************change
    # get the clip read for the breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        r_read = [
            read for read in bamf.fetch(r[0], r[1], r[1]+1)  # chr  coord  coord+1
            if read.reference_end == r[1]+1
            and check_read_clip(read, 'soft_right')]

    # get the sequence bases within the breakpoint
    start =  r[1] - seq_len + 1
    if  start < 0:
        start = 0
    r_seq = get_consensus_sequence(bamfile, r[0], start, r[1]+1)

    # get the clipped reads of the breakpoint
    r_clipped = [read.query_sequence[read.query_alignment_end:]
                 for read in r_read]
    r_clipped_name = [read.query_name for read in r_read]

    return [r_seq, r_clipped, r_clipped_name]

#------------------------------------------------------------------------------#
def L2R(l_clipped, r_seq):
    # compare b_clipped to a_seq (L2R)
    _clipped = [ax[::-1] for ax in l_clipped]
    _seq = r_seq[::-1]
    return (_clipped, _seq)

def R2L(r_clipped, l_seq):
    # compare b_clippled to a_seq (R2L)
    ######################## no necessary？？？############################
    _clipped = [bx for bx in r_clipped]
    _seq = l_seq
    return (_clipped, _seq)

def L2L(l_clipped, l_seq):
    # compare b_clippled to a_seq (L2L)
    _clipped = [rev_compl(bx) for bx in l_clipped]
    _seq = l_seq
    return (_clipped, _seq)

def R2R(r_clipped, r_seq):
    # compare b_clippled to a_seq (R2R)
    _clipped = [rev_compl(bx)[::-1] for bx in r_clipped]
    _seq = r_seq[::-1]
    return (_clipped, _seq)

#------------------------------------------------------------------------------#
def perfect_match(_clipped, _seq, num = 0, offset = 0, min_nt = 2):
    """
    perfect match between clipped reads and breakend sequence
    for whole length

    min_nt : minimal nt match is better to be larger than 2,
        1 nt match is 25% change by random,
        2 nt match is 0.0625 change by random,
        3 nt match is 0.0156 change by random,
    """
    if (offset > 0):
        m = [re.search(_x[offset:], _seq)
             for _x in _clipped if (len(_x) > offset) and (len(_x) >= min_nt)]
    else:
        m = [re.search(_x, _seq) for _x in _clipped if (len(_x) >= min_nt)]

    mm = [i.start() for i in m if i is not None]
    if (mm):
        if num > 0 :
            # get the most common split point
            match_start, scount =  Counter(mm).most_common()[0]
            scount += num
            return (match_start, scount)
        else:
            return Counter(mm).most_common()[0]
    else:
        return (0, 0)

#------------------------------------------------------------------------------#
def partial_match(_clipped, _seq, num = 0, min_nt = 2):
    """
    perfect match between clipped reads and breakend sequence
    for last n bases, mean the leading bases can be absent for potential
    insertions between two breakend

    """
    m = [SequenceMatcher(None, _x, _seq, autojunk = False)\
        .find_longest_match(0, len(_x), 0, len(_seq))
        for _x in _clipped if (len(_x) >= min_nt)]
    _clipped_len = [len(_x) for _x in _clipped]
    mm = [m[i].b for i in range(len(m)) \
        if (m[i].size + m[i].a == _clipped_len[i])]
    if (mm):
        if num > 0:
            # get the most common split point
            match_start, scount =  Counter(mm).most_common()[0]
            scount += num
            return (match_start, scount)
        else:
            return Counter(mm).most_common()[0]
    else:
        return (0, 0)

#------------------------------------------------------------------------------#
def fuzzy_match(_clipped, _seq, num = 0, offset = 0, min_nt = 2, error_rate = 0.1):
    """fuzzy match between clipped reads and breakend sequence
    use levenshtein distance to search complete sequence in _clipped
    that match part of _seq from start with error rate of 10%

    """
    match_start = 0
    if (offset > 0): # offset > 0 no overlap
        m = [levenshtein_distance(_x[offset:], _seq)
             for _x in _clipped if (len(_x) > offset) and (len(_x) >= min_nt)]
    else:
        match_start = -offset
        _sub_seq = _seq[match_start:]
        m = [levenshtein_distance(_x, _sub_seq)
             for _x in _clipped if (len(_x) >= min_nt)]
    if num > 0:
        return (match_start, sum(np.array(m) <= error_rate) + num)
    else:
        return (match_start, sum(np.array(m) <= error_rate))

#------------------------------------------------------------------------------#
def levenshtein_distance(str1, str2, rate = True, _sub = 1, _del = 1, _ins = 1):
    """levenshtein distance
    calculate the distance between two strings
    """
    # init zero
    rows = len(str1) + 1
    cols = len(str2) + 1
    dist = np.zeros((rows, cols), dtype = int)

    # init first row and column
    for i in range(1, rows):
        for j in range(1, cols):
            dist[i][0] = i
            dist[0][j] = j

    # dynamic programming
    for col in range(1, cols):
        for row in range(1, rows):
            if str1[row - 1] == str2[col - 1]:
                mm = 0 # match
            else:
                mm = _sub # substitution
            dist[row][col] = min( # minimal penalty
                dist[row - 1][col - 1] + mm, # substitution
                  dist[row - 1][col] + _del, # deletion
                dist[row][col - 1] + _ins, # insertion
            )
    if (rate):
        return dist[row].min()/len(str1)
    else:
        return dist[row].min()

def combinations_bp(bp_cand_df):
    '''
    combination cand bp pair from bp cand.
    '''
    bpp_cand_list = []
    for i in set(bp_cand_df['Group']):
        df = bp_cand_df[bp_cand_df['Group'] == i].sort_values(['Chrom', 'Coord', 'Clip'])
        if i != 2:
            bp_cand_list = [x[1:4] for x in df.itertuples()]
            for i,j in combinations(range(len(bp_cand_list)),2):
                bpp_ = (bp_cand_list[i], bp_cand_list[j]) 
                if bpp_ not in bpp_cand_list:
                    bpp_cand_list.append(bpp_)
        elif i == 2:
            for chrom in set(df['Chrom']):
                bp_cand_list = [x[1:4] for x in df[df['Chrom'] == chrom]\
                                .sort_values(['Chrom', 'Coord', 'Clip']).itertuples()]
                for i,j in combinations(range(len(bp_cand_list)),2):
                    bpp_ = (bp_cand_list[i],  bp_cand_list[j]) 
                    if bpp_ not in bpp_cand_list:
                        bpp_cand_list.append(bpp_)
    return bpp_cand_list
    

#------------------------------------------------------------------------------#
# main function to search for breakpoint duo in a large list of breakpoints
#------------------------------------------------------------------------------#
def get_breakpoint_duo(bamfile, bp_cand_df, \
    seq_len = 50, seed_len = None, min_nt = 5, rm_bpp = []):
    """
    get breakpoint duo from a list of breakpoints
    use pairwise assembly

    parameters:
    seq_len - # bases within breakend
    seed_len - # of bases up and down stream of the breakend
        in assembled sequence. up: seed_len; down: seed_len respectively
    """
   
    # get all bp cand
    bplist = [x[1:4] for x in bp_cand_df.drop_duplicates(['Chrom', 'Coord', 'Clip']).itertuples()]
  
    if len(set(bp_cand_df['Group'])) > 1:
        # find the breakpoint reads sequence
        bpass = {}
        for a in bplist:
            bpass[a] = breakpoint_assemble(bamfile, a, seq_len)

        # (a, b)
        bpp_cand_list = combinations_bp(bp_cand_df)  # may be not needed

        if rm_bpp:
            bpp_cand_list =  [i for i in bpp_cand_list if i not in rm_bpp]

        op = []

        for i, j in bpp_cand_list:
            mtf, offset, ab = breakpoint_matcher(bpass[i], bpass[i], \
                i[2] + j[2], seq_len, seed_len, min_nt)
            if (mtf):
                op.append(i + j + (mtf, offset, ab)) 
    else:
        # find the breakpoint reads sequence
        bpass = [breakpoint_assemble(bamfile, a, seq_len) for a in bplist]

        op = []
        for i, j in combinations(range(len(bplist)), 2):
            if (bplist[i], bplist[j]) not in rm_bpp:
                mtf, offset, ab = breakpoint_matcher(bpass[i], bpass[j], \
                    bplist[i][2] + bplist[j][2], seq_len, seed_len, min_nt)
                if (mtf):
                    op.append(bplist[i] + bplist[j] + (mtf, offset, ab))
    return op

def breakpoint_assemble(bamfile, a, seq_len = 50):
    """
    assemble sequence arround clipped breakpoint
    """
    if (a[2] == 'L'):
        a_seq, a_clipped, a_clipped_name = left_clip(bamfile, a, seq_len)
        if (a_clipped != []):
            # a_asmb = ''
            a_asmb = consensus_from_listofseq(a_clipped, a[2]) + a_seq
        else:
            a_asmb = ''
    else:
        a_seq, a_clipped, a_clipped_name = right_clip(bamfile, a, seq_len)
        if (a_clipped != []):
            # a_asmb = ''
            a_asmb = a_seq + consensus_from_listofseq(a_clipped, a[2])
        else:
            a_asmb = ''
    return a_asmb

#------------------------------------------------------------------------------#
def consensus_from_listofseq(rc, pattern = 'R'):
    """
    generate consensus sequence from a list of sequences
    """
    op = []
    length = max([len(x) for x in rc])

    if (pattern == 'R'):
        for j in range(length):
            op.append(Counter([s[j]
                      for s in rc if len(s) > j]).most_common(1)[0][0])
        return "".join(op)
    else:
        for j in range(length):
            op.append(Counter([s[::-1][j]
                      for s in rc if len(s) > j]).most_common(1)[0][0])
        return "".join(op)[::-1]


#------------------------------------------------------------------------------#
def breakpoint_matcher(a_ass, b_ass, pattern,
                       seq_len = 50, seed_len = None, min_nt = 5):
    """
    breakpoint match given two sequence

    output:
        offset:
            0  - perfect join ['ATTGC', 'GTCAT'] -> 'ATTGCGTCAT'
            1  - 1 bp gap (insertion) ['ATTGC', 'GTCAT'] -> 'ATTGCAGTCAT'
            -1 - 1 bp overlap ['ATTG', 'GTCAT'] -> 'ATTGTCAT'
    """
    # format assembled sequences to match
    # b->a = 5'->3'pattern
    if (pattern == 'RL'):
        a_ass, b_ass = b_ass, a_ass
    elif (pattern == 'RR'):
        a_ass = rev_compl(a_ass)
    elif (pattern == 'LL'):
        b_ass = rev_compl(b_ass)

    if (seed_len is not None): # match seed sequence first to increase speed
        seed_seq = b_ass[seq_len-seed_len:seq_len+seed_len]
        # seed_match = re.search(seed_seq, a_ass)
        # if (seed_match is not None):
        if (seed_seq in a_ass):
            m = SequenceMatcher(None, a_ass, b_ass, autojunk = False)
            m = m.find_longest_match(0, len(a_ass), 0, len(b_ass))  # longest continuous match
            ab = b_ass[0: m.b].lower() \
                + b_ass[m.b: m.b + m.size] \
                + a_ass[m.a + m.size:].lower() # sequence of the joint point
            is_paired = seq_len - min_nt >= m.b \
                and seq_len + min_nt <= m.b + m.size - 1
            offset = (len(a_ass) - seq_len - m.a) - (seq_len - 1 - m.b) - 1
        else:
            is_paired = False
            offset = None
            ab = None
    else: # don't use seed
        m = SequenceMatcher(None, a_ass, b_ass, autojunk = False)
        m = m.find_longest_match(0, len(a_ass), 0, len(b_ass))
        ab = b_ass[0: m.b].lower() \
            + b_ass[m.b: m.b + m.size] \
            + a_ass[m.a + m.size:].lower() # sequence of the joint point
        is_paired = seq_len - min_nt >= m.b \
            and seq_len + min_nt <= m.b + m.size - 1
        offset = (len(a_ass) - seq_len - m.a) - (seq_len - 1 -m.b) - 1

    return (is_paired, offset, ab)

#------------------------------------------------------------------------------#

def remove_duplicates(bp_duo, dist = 100, col = 1):
    '''
    merge the two bp pairs with the same start when two have the closer end.
    bp_duo: data frame of all bp pairs
    dist: merge only the distance of two ends less than dist
    '''
    bp_duo_new = []
    if col == 1:
        groupby = ['Chrom1','Coord1', 'Clip1']
        ind = 'Coord2'
    else:
        groupby = ['Chrom2','Coord2', 'Clip2']
        ind = 'Coord1'
    for i in bp_duo.groupby(groupby):
        df = i[1].sort_values(ind)
        # only do when has nearly end
        if df.shape[0] > 1 and (np.diff(df[ind]) < dist).any():
            final_result = []
            # cumulative_count = 0
            df = df.sort_values(by = ind)
 
            dist1 = dist
            for index, row in df.reset_index(drop=True).iterrows():
                if index == 0:
                    # cumulative_count = row['Count']
                    prev_row = row
                    t = 0
                else:
                    distance = row[ind] - prev_row[ind]
                    if (distance > dist) and  (distance > dist1 ) and (t != 0):
                        t = 0
                        
                    if distance <= dist or (distance <= dist1 and dist1 != dist):
                        t += 1
                        # keep the row with largest count 
                        if prev_row['Count'] >= row['Count']:   
                            dist1 = dist + distance
                        else:
                            prev_row = row
                            dist1 = dist
                    else:
                        final_result.append(prev_row)
                        prev_row = row
                        dist1 = dist
            final_result.append(prev_row)
            bp_duo_new.append(pd.DataFrame(final_result))
        else:
            bp_duo_new.append(i[1])

    return pd.concat(bp_duo_new)



#------------------------------------------------------------------------------#
def get_bppair(bamfile, bp_cand_df, bp_duo_bam,\
    seq_len = 50, seed_len = 5, min_nt = 5,
    match_method = 'fuzzy_match', multi_process = False, num_cores = 8):
    """
    get the bppairs from bp_cand_stats (a list of bps)

    parameters:
    seq_len - # bases within breakend
    seed_len - # of bases up and down stream of the breakend
        in assembled b sequence. up: seed_len; down: seed_len respectively
    """
    
    # # example of element of bplist： ('chr7', 54817063, 'L')
    # bplist = [x[1:4] for x in bp_cand_df_sorted.itertuples()]
    if bp_duo_bam.shape[0] > 10000:
        multi_process = True
        num_cores = 20

    if not bp_duo_bam.empty:
        rm_bpp = bp_duo_bam.index.tolist()  #  c
        # rm_bpp = [(i[1:4],i[4:7]) for i in bp_duo_bam.itertuples()]
    else:
        rm_bpp = []
        
    # get the breakpoint pairs
    # note: searching use a different (shorter) seq_len parameter
    # to increase running efficiency
    if not bp_cand_df.empty:  # this is can remove
        #### bpp_cand_list is the bp_cand_df
        bp_cand_df_sorted = bp_cand_df.sort_values(['Chrom', 'Coord', 'Clip'])
        bpduo = get_breakpoint_duo(bamfile, bp_cand_df_sorted, seq_len, seed_len, min_nt, rm_bpp = rm_bpp)  # consuming time 
    else:
        bpduo = []

    # count the # of supporting reads (clip reads)
    # note: counting use a different (longer) seq_len parameter
    # to ensure the reads are fully coverred in the breakend
    def run_join_breakpoint(bamfile, match_method, type, t1, row):
        if type == 1:
            # logg.info(f'bpduo: {row[0:6]} of {len(t1)}')
            bp_count = tuple(join_breakpoint(bamfile, row[0:3], row[3:6], rm_clipped_name = None, \
                    offset = row[7], match_method = match_method))
            if bp_count[0][1] > 0 and bp_count[1][1] > 0 and bp_count[0][0] == bp_count[1][0]:
                sum_count = bp_count[0][1] + bp_count[1][1]
                t1.append(list(row[0:6] + (sum_count,) + row[7:]))
            return t1
        else: 
            # logg.info(f'bpduobam: {row[1:7]} of {len(t1)}')
       
            bp_count = join_breakpoint(bamfile, row[1:4], row[4:7], rm_clipped_name = row[-1], \
                    offset = row[8], match_method = match_method)
            if bp_count[0][0] == bp_count[1][0]:
                sum_count = bp_count[0][1] + bp_count[1][1]
                if sum_count < row[7]:  # raw count bp_pair
                    sum_count = row[7]
                t1.append(list(row[1:4] + row[4:7] + (sum_count,) + row[8:10]))   
            return t1
    
    if multi_process:
        t1 = []
        if not bp_duo_bam.empty:
            pool = Pool(num_cores)
            func = partial(run_join_breakpoint, bamfile, match_method, 2, t1)
            pool.map(func, bp_duo_bam.itertuples())  # [bp_duo_bam['Count']>1]
        if bpduo:
            pool = Pool(num_cores)
            func2 = partial(run_join_breakpoint, bamfile, match_method, 1, t1)
            pool.map(func2, bpduo)
    else:
        t1 = []
        if bpduo:
            tmp = 0
            for row in bpduo:  # [:10]
                logg.info(f'bpduo: {tmp} of {len(bpduo)}')
                t1 = run_join_breakpoint(bamfile, match_method, type = 1, t1 = t1, row = row)
                tmp += 1
        if not bp_duo_bam.empty:
            tmp = 0
            for row in bp_duo_bam[bp_duo_bam['Count'] >= 1].itertuples():  # >2 改为1
            # for row in bp_duo_bam_1[bp_duo_bam_1['Count'] >= 1].itertuples():  # >2 改为1
                logg.info(f'bpduobam: {tmp} of {bp_duo_bam.shape[0]}')
                t1 = run_join_breakpoint(bamfile, match_method, type = 2, t1 = t1, row = row)
                tmp += 1

    colnames = ["Chrom1", "Coord1", "Clip1",
                "Chrom2", "Coord2", "Clip2", 'Count', 'offset', 'Seq']

    bp_pair_df = pd.DataFrame(t1, columns = colnames).reset_index(drop=True)
    if bp_pair_df.shape[0] > 1:
        bp_pair_df = remove_duplicates(bp_pair_df, dist = 2, col = 1)
        bp_pair_df = remove_duplicates(bp_pair_df, dist = 2, col = 2)
    return bp_pair_df.sort_values(
        'Count', ascending = False).reset_index(drop=True)

#------------------------------------------------------------------------------#
# extract break joint sequences
#------------------------------------------------------------------------------#
def left_clip_seq(bamfile, l, seq_len = 100):
    # get the clip read for breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        l_read = [
            read for read in bamf.fetch(l[0], l[1], l[1]+1) \
            if check_read_clip(read, 'soft_left', seq = True) \
            and read.reference_start == l[1]]

    # get the clipped reads of the breakpoint
    l_clipped_max = max([read.query_alignment_start for read in l_read])
    l_clipped = [' ' * (l_clipped_max-read.query_alignment_start) +
        read.query_sequence for read in l_read]
    l_clipped.sort()

    # get the sequence bases within the breakpoint
    l_seq = get_consensus_sequence(bamfile, l[0], l[1], l[1]+seq_len)
    l_seq = ' ' * (l_clipped_max) + l_seq

    return [l_seq] + l_clipped, l_clipped_max

#------------------------------------------------------------------------------#
def right_clip_seq(bamfile, r):
    # get the clip read for the breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        r_read = [
            read for read in bamf.fetch(r[0], r[1], r[1]+1) \
            if check_read_clip(read, 'soft_right', seq = True) \
            and read.reference_end == r[1]+1]

    # get the clipped reads of the breakpoint
    r_clipped_max = max([read.query_alignment_end for read in r_read])
    r_clipped = [' ' * (r_clipped_max-read.query_alignment_end) +
        read.query_sequence for read in r_read]
    r_clipped.sort(reverse = True)

    # get the sequence bases within the breakpoint
    r_seq = get_consensus_sequence(bamfile, r[0], r[1]-r_clipped_max+1, r[1]+1)

    return [r_seq] + r_clipped, r_clipped_max

#------------------------------------------------------------------------------#
def get_alignment(bamfile, a, b, offset = 0):
    """get alignment of bp pairs
    """
    if (a[2] == 'R'):
        a_seq, a_clip_max = right_clip_seq(bamfile, a)
    else:
        a_seq, a_clip_max = left_clip_seq(bamfile, a)

    if (b[2] == 'R'):
        b_seq, b_clip_max = right_clip_seq(bamfile, b)
    else:
        b_seq, b_clip_max = left_clip_seq(bamfile, b)

    pattern = a[2] + b[2]

    if (pattern == 'RL'):
        b_seq = [' ' * (a_clip_max - b_clip_max + offset) + read
                 for read in b_seq]
    elif (pattern == 'LR'):
        a_seq, b_seq = b_seq, a_seq
        a_clip_max, b_clip_max = b_clip_max, a_clip_max
        b_seq = [' ' * (a_clip_max - b_clip_max + offset) + read
                 for read in b_seq]
    elif (pattern == 'LL'):
        a_clipped_max = max([len(z) for z in a_seq])
        a_seq = [' ' * (a_clipped_max - len(z)) + rev_compl(z.strip())
                 for z in a_seq]
        b_seq = [' ' * (a_clipped_max - a_clip_max - b_clip_max + offset)
                 + read for read in b_seq]
    else:
        b_clipped_max = max([len(z) for z in b_seq])
        b_seq = [' ' * (b_clipped_max - len(z)) + rev_compl(z.strip())
                 for z in b_seq]
        b_seq = [' ' * (len(a_seq[0]) + offset - b_seq[0].count(' '))
                 + read for read in b_seq]
    # Correct number of space    
    tmp = len(b_seq[0]) - len(b_seq[0].strip())  - (len(a_seq[0]) + offset) 
    if tmp > 0:
        b_seq = [bseq[tmp:] for bseq in b_seq]
    return [a_seq[0]] + [b_seq[0]] + [] + a_seq[1:] + b_seq[1:]

#------------------------------------------------------------------------------#
def ouput_alignment(bp_pair, bamfile, output_align_dir, circ_anno):
    """output alignment on pairs of breakpoints
    """
    if not circ_anno.empty:
        AmpliconIDs = circ_anno['AmpliconID'].unique()
        for i in AmpliconIDs:
            output_align = output_align_dir + '/' + i + '.tsv'
            circ_anno_sub = circ_anno[circ_anno['AmpliconID'] == i]
            # get the bp_pair within  circ_anno
            if not os.path.exists(output_align_dir):
                os.makedirs(output_align_dir)
            ## get breakpoints 
            circ_anno_sub_copy = circ_anno_sub.copy()
            circ_anno_sub_copy.loc[circ_anno_sub_copy['Strand'] == '-', ['Start', 'End']] = circ_anno_sub_copy.loc[circ_anno_sub_copy['Strand'] == '-', ['End', 'Start']].values
            with open(output_align, 'w') as f:
                for id in range(len(circ_anno_sub_copy)):
                    current_row = circ_anno_sub_copy.iloc[id]
                    next_row = circ_anno_sub_copy.iloc[(id + 1) % len(circ_anno_sub_copy)]  #
                    if current_row[2] < next_row[1]:
                        index = ((bp_pair['Coord1'] == current_row[2]) & (bp_pair['Chrom1'] == current_row[0])) & \
                            ((bp_pair['Coord2'] == next_row[1]) & (bp_pair['Chrom2'] == next_row[0]))
                    else:
                        index = ((bp_pair['Coord2'] == current_row[2]) & (bp_pair['Chrom2'] == current_row[0])) & \
                            ((bp_pair['Coord1'] == next_row[1]) & (bp_pair['Chrom1'] == next_row[0]))
                    bp_pair_sub = bp_pair[index]
                
                    for row in bp_pair_sub.itertuples():
                        if (row[9] != 'PE_Support'):
                            a = row[1:4]
                            b = row[4:7]
                            offset = row[8]
                            output_align_tmp_file = output_align_dir + '/' +  a[0] + str(a[1]) + a[2] + b[0] +  str(b[1]) + b[2] +'_tmp.tsv'
                            if not os.path.exists(output_align_tmp_file):
                                try:
                                    alignment = get_alignment(bamfile, a, b, offset)
                                    with open(output_align_tmp_file, 'w') as fin:
                                        fin.write('>' + '\t'.join(map(str, list(row[1:9]))) + "\n")
                                        for read in alignment:
                                            fin.write("%s\n" % read)
                                        fin.write("\n\n")
                                    with open(output_align_tmp_file, 'r') as fin:
                                        f.write(fin.read())
                                except:
                                    print(f'warings: {row}')
                            else:
                                with open(output_align_tmp_file, 'r') as source_file:
                                    f.write(source_file.read())
        cmd = 'rm -rf ' + output_align_dir + '/' + '*' + '_tmp.tsv'
        os.system(cmd)         

        

#------------------------------------------------------------------------------#
# search for breakpoint pairs candidates from improper paired reads
#------------------------------------------------------------------------------#
def get_bppair_peread(bamfile, bp_cand_df, \
    seq_len = 500, min_reads = 5):
    """
    get the bppairs from improper paired reads

    parameters:
    seq_len : # bases within breakend to search for
        improper paired reads
    min_reads : # paired reads supporting the breakpoint
    """

    # get the improper paired reads for each breakend
    bplist = []
    for row in bp_cand_df.itertuples():
        if (row[3] == 'L'):
            a = [row[1], row[2], row[2] + seq_len]
        elif (row[3] == 'R'):
            a = [row[1], row[2] - seq_len, row[2]]
        with pysam.AlignmentFile(bamfile, "rb") as bamf:
            a_read = [
                read.query_name for read in bamf.fetch(a[0], a[1], a[2])
                if (check_read(read) and (read.is_proper_pair == False))]
        if (len(a_read) >= min_reads):
            bplist.append([row[1], row[2], row[3], a_read])

    # search pairs
    op = []
    for i,j in combinations(range(len(bplist)), 2):
        if (bplist[i][0] == bplist[j][0] and
            abs(bplist[i][1] - bplist[j][1]) < seq_len):
            isc = 0
        else:
            isc = len(misc.intersect(bplist[i][3], bplist[j][3]))

        if (isc >= min_reads):
            op.append(bplist[i][0:3] + bplist[j][0:3] + [isc, 0, 'PE_Support'])
    colnames = ["Chrom1", "Coord1", "Clip1",
                "Chrom2", "Coord2", "Clip2", 'Count', 'offset', 'Seq']
    bp_pair_df = pd.DataFrame(op, columns = colnames)

    return bp_pair_df.sort_values(
        'Count', ascending=False).reset_index(drop=True)


#------------------------------------------------------------------------------#
def bp_cn_boundary(cn_amp, bp_cand_stats = None):
    """
    get breakpoints from breakpoint candidates and amplified segments
    """
    if (bp_cand_stats is not None):
        ## make breakpoints dataframe
        bplist = [i[1:4] for i in bp_cand_stats.itertuples()]
        op = []
        for row in cn_amp.itertuples():
            if ((row[1], row[2], 'L') not in bplist):
                op.append([row[1], row[2], 'L', True, row[5], 0, 0])
            if ((row[1], row[3], 'R') not in bplist):
                op.append([row[1], row[3], 'R', True, row[6], 0, 0])
        cn_seg_df = pd.DataFrame(op, columns = bp_cand_stats.columns)

        # merge and output
        df = pd.concat([bp_cand_stats, cn_seg_df])
        df = df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
        return df.reset_index(drop = True)
    else:
        op = []
        for row in cn_amp.itertuples():
            op.append([row[1], row[2], 'L', True, row[5], 0, 0])
            op.append([row[1], row[3], 'R', True, row[6], 0, 0])
        columns = ['Chrom', 'Coord', 'Clip', 'CleanBP',
                   'ClipDepth', 'InDepth', 'OutDepth']
        cn_seg_df = pd.DataFrame(op, columns = columns)
        df = cn_seg_df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
        return df.reset_index(drop = True)