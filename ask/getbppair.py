import os
import gzip
import pysam
import pickle
import pandas as pd
import sys
import os
from collections import Counter, defaultdict
import time
#------------------------------------------------------------------------------#


##################################################################
#  this code is used to find bp-pair from bamfile.
##################################################################
def check_read(read, mapq = 20, nmmax = 1):
    """
    true if the read meets certian conditions
    """
    # read.get_tag('NM') <= nmmax ？？？？？？？
    return not read.is_unmapped \
        and read.mapping_quality >= mapq \
        and read.get_tag('NM') <= nmmax \
        and not read.is_duplicate

def rev_compl(seq):
    nn = defaultdict(
        lambda: 'N', {
            'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
            'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'})
    return ''.join([nn[a] for a in seq][::-1])

def get_clip_read(read):
    '''
    get clip read
    '''
    strand = '+'
    if read.is_reverse:
       strand = '-'
    if (read.cigartuples[0][0] == 4): # soft clip reads
        return (read.reference_name, \
                        read.reference_start, strand , "L", read.cigarstring,\
                            [read.is_unmapped,read.mapping_quality,read.get_tag('NM'),read.is_duplicate])
    elif (read.cigartuples[0][0] == 5): # hard clip reads
        return (read.reference_name, \
                        read.reference_start,strand, "L", read.cigarstring,\
                            [read.is_unmapped,read.mapping_quality,read.get_tag('NM'),read.is_duplicate])
    # right clip reads
    if (read.cigartuples[-1][0] == 4): # soft clip reads
        return (read.reference_name, \
                        read.reference_end-1, strand, "R", read.cigarstring,\
                            [read.is_unmapped,read.mapping_quality,read.get_tag('NM'),read.is_duplicate])
    elif (read.cigartuples[-1][0] == 5): # hard clip reads
        return (read.reference_name, \
                        read.reference_end-1, strand, "R", read.cigarstring,\
                            [read.is_unmapped,read.mapping_quality,read.get_tag('NM'),read.is_duplicate])


def join_breakpoint_(read, read_s, suppl, a, b ):
    diff_strand = a[2] != b[2]
    if a[3] + b[3] == 'LL':
        if b[0] in suppl[0] and b[2] == suppl[2]:  # and a[1] == int(suppl[1])-1 
            if not diff_strand:
                if read_s.cigartuples[0][1] > read.cigartuples[0][1]: 
                    # this refers a is R clip, a-->b
                    if read.cigartuples[-1][0] == 4:
                        a = (a[0], a[1] + read.cigartuples[1][1] - 1, a[2], 'R', a[4], a[5])
                    which_trans = 'a'
                elif read_s.cigartuples[0][1] < read.cigartuples[0][1]:
                    if read_s.cigartuples[-1][0] == 5:
                        b = (b[0], b[1] + read_s.cigartuples[1][1] - 1, b[2], 'R', b[4], b[5])
                    which_trans = 'b' 
            else:
                which_trans = None
            bp_pair = a + b
            ba, offset, b_seq, a_seq = L2L_seq(read, read_s, diff_strand, which_trans)
    elif a[3] + b[3] == 'RL':
        if b[0] in suppl[0] and b[2] == suppl[2]: # and b[1] == int(suppl[1])-1 
            if diff_strand:
                if read_s.cigartuples[-1][0]==5:  
                    b = (b[0], b[1] + read_s.cigartuples[1][1] - 1, b[2], 'R', b[4], b[5])
                    which_trans = 'b'
            else:
                which_trans = None
            bp_pair = a + b
            ba, offset, b_seq, a_seq = R2L_seq(read, read_s, diff_strand, which_trans)
    elif a[3] + b[3] == 'RR':
        if b[0] in suppl[0] and b[2] == suppl[2]:
            bp_pair = a + b
            if diff_strand:
                ba, offset, b_seq, a_seq = R2R_seq(read, read_s, diff_strand)
    elif a[3] + b[3] == 'LR':
        if b[0] in suppl[0] and b[2] == suppl[2]:
            if diff_strand:
                if read.cigartuples[-1][0] == 4:
                    a = (a[0], a[1] + read_s.cigartuples[1][1] - 1, a[2], 'R', a[4], a[5])
                    which_trans = 'a'
            else:
                which_trans = None
            bp_pair = a + b
            ba, offset, b_seq, a_seq = L2R_seq(read, read_s, diff_strand, which_trans)
    return a, b, bp_pair, ba, offset, b_seq, a_seq
#-----------------------------------------------------------------
#  get bp-pair sequence 
#-----------------------------------------------------------------
def L2R_seq(read, read_s, diff_strand = False, which_trans = None):
    # check wheather query_sequence is euqal with query_alignment_sequence (actually euqal)
    # s1 = read.cigartuples[0][1]
    # e1 = read.cigartuples[1][1]
    # record_query = read.query_sequence[s1:s1+e1]
    # is_euqal = record_query == read.query_alignment_sequence  # can remove
    if diff_strand:  # r:+, r_s: -
        # if is_euqal: # can remove
        #     b_seq = read.query_alignment_sequence  
        # else:
        #     b_seq = record_query
        #     a_seq = rev_compl(read_s.query_alignment_sequence)
        if which_trans == 'a':
            a_seq = rev_compl(read_s.query_alignment_sequence)
            b_seq = read.query_alignment_sequence
            ## mapped length of read1 - clipped length of supplementary (>0: overlap)
            offset = read.cigartuples[0][1] +  read.query_alignment_length  - read_s.cigartuples[-1][1]   # if b is L to R, s_offset=1
        else:
            print(f'{read.query_name}: this is impossible!!!')
            offset = 3
            a_seq = None
            b_seq = None
    else:
        a_seq = read.query_alignment_sequence
        b_seq = read_s.query_alignment_sequence
        # if 'I' in read_s.cigarstring:
        offset = read_s.query_alignment_length - read.cigartuples[0][1]
    if offset <= 0: 
        ba = b_seq + '-' * (-offset) + a_seq 
    elif offset > 0: 
        ba = b_seq + a_seq[offset:]
    return ba, offset, b_seq, a_seq

def L2L_seq(read, read_s, diff_strand = False, which_trans = None):
    ### opt 2
    ## mapped length of supplementary - clipped length of read1 ( >0: overlap )
    # check wheather query_sequence is euqal with query_alignment_sequence (actually euqal)
    a_seq = read.query_alignment_sequence
    b_seq = read_s.query_alignment_sequence
    if diff_strand:
        # #2: total hard clip bp - left clip
        offset = read_s.get_cigar_stats()[0][5] - read_s.cigartuples[0][1]\
              + read_s.query_alignment_length - read.cigartuples[0][1]
        b_seq = rev_compl(read_s.query_alignment_sequence)
    else:
        if which_trans == 'b': # 
            offset = read_s.cigartuples[0][1] + read_s.query_alignment_length - read.cigartuples[0][1]
        elif which_trans == 'a': # read: 4S45M27S, reads: 30H41M5H
            offset = read.query_alignment_length + read.cigartuples[0][1] - read_s.cigartuples[0][1]
            b_seq, a_seq = a_seq, b_seq
    if offset > 0:
        ba = b_seq  + a_seq[offset:] 
    elif offset <= 0:
        ba = b_seq + '-' * (-offset) + a_seq                    
    return ba, offset, b_seq, a_seq

def R2L_seq(read, read_s, diff_strand = False, which_trans = None):
    if diff_strand:  # r:+, r_s: -
        if which_trans == 'b':
            a_seq = rev_compl(read_s.query_alignment_sequence)
            b_seq = read.query_alignment_sequence
            ## mapped length of read1 - clipped length of supplementary (>0: overlap)
            offset = read.query_alignment_length - read_s.cigartuples[-1][1] 

        else: 
            print(f'{read.query_name}: this is impossible!!!')
    else:
        a_seq = read_s.query_alignment_sequence
        b_seq = read.query_alignment_sequence
        offset = read.query_alignment_length - read_s.cigartuples[0][1]  #  len(b_seq): for '7M4I32M 31S'
        # actuallly ba = ab
    if offset > 0:
        ba = b_seq + a_seq[offset:]
    elif offset <= 0: # check
        ba = b_seq + '-' * (-offset) + a_seq 
    return ba, offset, b_seq, a_seq

def R2R_seq(read, read_s, diff_strand = False): 
    if  diff_strand:
        offset = read.query_alignment_length - read_s.cigartuples[-1][1] #read_s.cigartuples[1][1]
        b_seq = read.query_alignment_sequence 
        a_seq = rev_compl(read_s.query_alignment_sequence)
        if offset <= 0:
            ba = b_seq + '-' * (-offset) + a_seq
        elif offset > 0:
            ba = b_seq  + a_seq[offset:]
    else: # check
        print(f'{read.query_name}: this is impossible!!!')
        offset = None
        if read.cigartuples[0][1] >= read_s.cigartuples[0][1]:
            b_seq = read.query_alignment_sequence
            a_seq = ' '
            ba = read.query_alignment_sequence
        else:
            b_seq = rev_compl(read_s.query_alignment_sequence)
            a_seq = ''
            ba = b_seq
    return ba, offset, b_seq, a_seq



#-----------------------------------------------------------------
#  get reads with SA
#-----------------------------------------------------------------

#---------------------------------------no use---------------------------------------
def get_reads_with_SA(bamfile):
    '''
    get bp pair reads from bamfile
    '''
    reads_with_SA = {}
    with pysam.AlignmentFile(os.path.join(bamfile), "rb") as samfile:
  
        for read in samfile.fetch():
            if read.has_tag('SA') and check_read(read,mapq=0,nmmax=1):
                print(read.query_name)
             
                read_name = read.query_name
      
                if read_name not in reads_with_SA:
                    reads_with_SA[read_name] = {}
        
                # save the pair of read and supplementary read with SA
                if read.is_read1 and not read.is_supplementary:
                    if 'read1' not in reads_with_SA[read_name]:
                        reads_with_SA[read_name]['read1'] = []
                    reads_with_SA[read_name]['read1'].append(read)
                elif read.is_read1 and read.is_supplementary:
                    if 'read_s' not in reads_with_SA[read_name]:
                        reads_with_SA[read_name]['read1_s'] = []
                    reads_with_SA[read_name ]['read1_s'].append(read)
                elif read.is_read2 and not read.is_supplementary:
                    if 'read2' not in reads_with_SA[read_name]:
                        reads_with_SA[read_name]['read2'] = []
                    reads_with_SA[read_name]['read2'].append(read)
                elif read.is_read2 and read.is_supplementary:
                    if 'read2_s' not in reads_with_SA[read_name]:
                        reads_with_SA[read_name]['read2_s'] = []
                    reads_with_SA[read_name ]['read2_s'].append(read)
                else:
                    reads_with_SA[read_name]['single'].append(read)
    return reads_with_SA


def check_bp_pair(readp):
    if 'read1' in readp and 'read1_s' in readp:
        return 'read1'
    elif 'read2' in readp and 'read2_s' in readp:
        return 'read2'
    else:
        return None


def get_read_barcode(read, tags=("CB", "CR", "BX", "BC")):
    for tag in tags:
        if read.has_tag(tag):
            return read.get_tag(tag)
    return None


'''
bp_pair = 
('chrY', 58992343, '+', 'L', '17S58M', ['SRR8085207.430919', False, 0, 5, False],
 'chrY', 58991061, '+', 'R', '33M42H', ['SRR8085207.430919', False, 0, 1, False])
'''
bp_pair_list = []
def get_bp_pair(reads_with_SA, bp_pair_info_file = None, mapq = 20):
    bp_pair_dict = {}
    for readname, readp in reads_with_SA.items():
        # print(f'read: {readname},len is {len(readp)}.')
        if len(readp) == 2:
            key = check_bp_pair(readp)
            if key:
                read = readp[key][0]
                suppls = [j.split(',') for j in read.get_tag('SA').split(';') if j !='']  # only S
                for read_s in readp[key+'_s']:  # supplmentary read may be multiple.
                    if (read.mapping_quality >= mapq and read.get_tag('NM') <= 1) or (read_s.mapping_quality >= mapq and read_s.get_tag('NM') <= 1):
                        # read_s = readlen(readp) p[key+'_s'][0]
                        b = get_clip_read(read_s)
                        a = get_clip_read(read)
                        strand = '-'
                        if read_s.is_forward:
                            strand = '+'
                        index = [id for id, j in enumerate(suppls) if j[3] == b[4].replace('H','S') and j[2] == strand and j[0] == b[0]]
                        
                        try:
                            if bp_pair_info_file:
                                with open(bp_pair_info_file, 'a') as f:
                                    f.write(f"readname_{read.query_name}\tread_sname_{read_s.query_name}\t{a}\t{b}\tafter\t")  
                            a, b, bp_pair, ba, offset, b_seq, a_seq = join_breakpoint_(read,read_s, suppls[index[0]], a, b)
                                                        
                            #-----------
                            offset = -offset 
                            if bp_pair[1] < bp_pair[7]:
                                bpp = (bp_pair[:2]+ bp_pair[3:4], bp_pair[6:8] + bp_pair[9:10])
                            else:
                                # bp_pair[6:8] + bp_pair[9:10]+ bp_pair[:2]+ bp_pair[3:4]
                                bpp = (bp_pair[6:8] + bp_pair[9:10], bp_pair[:2]+ bp_pair[3:4])
                            # dkey = str(bpp)
                            readsbarcode = get_read_barcode(read)
                            if bpp not in bp_pair_dict:
                                bp_pair_dict[bpp] = list(bpp[0] + bpp[1])  + [1] + [offset] + [ba] + [[read.query_name]] + [[readsbarcode]]
                            else:
                                # if bpp in bp_pair_dict:
                                # if len(bp_pair_dict[bpp][-1].replace('-','')) < len(ba.replace('-','')):  # keep the longest reads
                                if bp_pair_dict[bpp][7] > offset:
                                    bp_pair_dict[bpp][8] = ba 
                                    bp_pair_dict[bpp][7] = offset 
                                bp_pair_dict[bpp][6] += 1       
                                bp_pair_dict[bpp][9].append(read.query_name)
                                bp_pair_dict[bpp][10].append(readsbarcode)
                            
                            #-----------

                            if bp_pair_info_file:
                                with open(bp_pair_info_file, 'a') as f:
                                    f.write(f"{readname}\t{a}\t{b}\t{offset}\t{read.query_sequence}\t{read.query_alignment_sequence}\t{read_s.query_alignment_sequence}\t{b_seq}\t{a_seq}\t{ba}\n")         
                        except:
                            print(f'{readname}: error!!!')
        # -------- if read1 and read2 both have supplemnetary------
        elif len(readp) == 4:
            key = check_bp_pair(readp)
            if key:
                for i in range(1,3):
                    read = readp['read' + str(i)][0]
                    suppls = [j.split(',') for j in read.get_tag('SA').split(';') if j !='']  # only S
                    for read_s in readp['read' + str(i) + '_s']:  # supplmentary read may be multiple.
                        if (read.mapping_quality >= mapq and read.get_tag('NM') <= 1) or (read_s.mapping_quality >= mapq and read_s.get_tag('NM') <= 1):    
                            a = get_clip_read(read)
                            b = get_clip_read(read_s)
                            # get the pair breakpoint 
                            strand = '-'
                            if read_s.is_forward:
                                strand = '+'
                            index = [id for id, j in enumerate(suppls) if j[3] == b[4].replace('H','S') and j[2] == strand and j[0] == b[0]]
                            try:
                                if bp_pair_info_file:
                                    with open(bp_pair_info_file, 'a') as f:
                                        f.write(f"readname_{read.query_name}\tread_sname{read_s.query_name}\t{a}\t{b}\tafter\t") 
                                a, b, bp_pair, ba, offset, b_seq, a_seq  = join_breakpoint_(read,read_s, suppls[index[0]], a, b)

                                #-----------
                                offset = -offset  # <0 represents overlap
                                if bp_pair[1] < bp_pair[7]:
                                    bpp = (bp_pair[:2]+ bp_pair[3:4], bp_pair[6:8] + bp_pair[9:10])
                                else:
                                    bpp = (bp_pair[6:8] + bp_pair[9:10], bp_pair[:2]+ bp_pair[3:4])
                                # dkey = str(bpp)
                                readsbarcode = get_read_barcode(read)
                                if bpp not in bp_pair_dict:
                                    bp_pair_dict[bpp] = list(bpp[0] + bpp[1]) + [1] + [offset] + [ba] + [[read.query_name]] + [[readsbarcode]]
                                else:
                                    # if bpp in bp_pair_dict:
                                    # if len(bp_pair_dict[bpp][-1].replace('-','')) < len(ba.replace('-','')):  # keep the longest reads
                                    if bp_pair_dict[bpp][7] > offset:
                                        bp_pair_dict[bpp][8] = ba 
                                        bp_pair_dict[bpp][7] = offset 
                                    bp_pair_dict[bpp][6] += 1       
                                    bp_pair_dict[bpp][9].append(read.query_name)
                                    bp_pair_dict[bpp][10].append(readsbarcode)
                                #-----------
                                # try2
                                # bp_pair_list.append(bp_pair[:2]+ bp_pair[3:4] + bp_pair[6:8] + bp_pair[9:10] + (offset,) + (ba,))
                                if bp_pair_info_file:
                                    with open(bp_pair_info_file, 'a') as f:
                                        f.write(f"{readname}\t{a}\t{b}\t{offset}\t{read.query_sequence}\t{read.query_alignment_sequence}\t{read_s.query_alignment_sequence}\t{b_seq}\t{a_seq}\t{ba}\n")         
                      
                            except:
                                print(f'{readname}: error!!!')
        # else:
        #     print(f'check: {readname}, len is {len(readp)}!!!')
    if len(bp_pair_dict) > 0:
        bp_pair_df_by_dict = pd.DataFrame(bp_pair_dict).T
        bp_pair_df_by_dict.index = list(bp_pair_dict.keys())
        columns = ['Chrom1','Coord1','Clip1','Chrom2','Coord2','Clip2','Count','offset','Seq', 'Readsname', 'Readsbarcode']
        bp_pair_df_by_dict.columns = columns
    else:
        bp_pair_df_by_dict = pd.DataFrame()
    return bp_pair_df_by_dict


