################################################################################
# APIs to AmpliconSeeK
################################################################################
import os
import sys
import pandas as pd
import numpy as np
import covbam
import pysam
import cndetect
import bpdetect
import bpjoint
import evidence
import bppair
import ggraph
import grange
import trackplot
import misc
import cbs
import getbppair
import time
import logger
import warnings

warnings.filterwarnings('ignore')

logg = logger.Logger()

import numpy as np
from memory_profiler import profile

#------------------------------------------------------------------------------#
def process_alignment(bamfile, gsfile = None, binsize = 10000,
                      mapq = 1, nmmax = 3,
                      mode = 'standard', sub_binsize = 1000,
                      seg_robost_quant = 0.5, SA_with_nm = True, downsample = 10):  # ******change
    """from bam file to amplified segments and breakpoint list

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    gsfile : genome size file (not yet in use)
    binsize : genomic bin size to count reads in
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)
    mode : bin read count mode
        'standard' : calculate total read counts in target bins
                     used for data with low bias, such as Input
        'bias' : calculate robust statistics in target bins
                 used for data with bias coverage, such as ChIP
    sub_binsize : smaller bin size within target bin size to
        calculate robust statistics
    seg_robost_quant : quantile as robust statistics
        if None, if None auto determine for different library size
    """
    filename = os.path.splitext(os.path.basename(bamfile))
    logg.info(f'process_alignment: {filename}')
    ## get bamfile
    bamfile_ = pysam.AlignmentFile(bamfile, 'r')
    bamfile_init = covbam.prebamfile(bamfile_, downsample = downsample)
    save_stat = os.path.splitext(bamfile)[0] + '_' + str(downsample) + '_covstat.tsv'
    if bamfile_init.bamfile_filesize < 50:
        bamfile_init.downsample_ratio = 1
    else:
        if os.path.exists(save_stat):
            stats = []
            with open(save_stat, 'r') as f:
                for row in f:
                    res = row.strip().split('\t')
                    for j in res:
                        stats.append(j)
            bamfile_init.basic_stats = tuple(stats[1:])
            bamfile_init.downsample_ratio = float(stats[0])
        else:
            bamfile_init.median_coverage(window_size = 10000)
            with open(save_stat, 'w') as f:
                f.write(f'{bamfile_init.downsample_ratio}\t')
                for i in bamfile_init.basic_stats:
                    f.write(f'{i}\t')
        bamfile_init.downsample_ratio = min(0.4, bamfile_init.downsample_ratio)
    # extract all clip reads
    if mode == 'standard':
        bp_all, reads_with_SA, bin_count = covbam.clip_from_bam(
            bamfile_init, mapq=mapq, nmmax=nmmax,
            SA_with_nm=SA_with_nm, count_binsize=binsize)
    else:
        bp_all, reads_with_SA = covbam.clip_from_bam(
            bamfile_init, mapq=mapq, nmmax=nmmax, SA_with_nm=SA_with_nm)
    if bamfile_init.downsample_ratio != 1:
        ds_bamfile_pathname = covbam.downsampled_bam_path(
            bamfile_init, bamfile_init.downsample_ratio)
        pysam.index(ds_bamfile_pathname)
    # output bedgraph
    clip_bg = covbam.clip2bedgraph(bp_all)

    # get read counts in genomic bins
    if mode != 'standard':
        bin_count = covbam.region_count_2pass(bamfile_init, gsfile,
            binsize = binsize, mapq = mapq, nmmax = nmmax,
            sub_binsize = sub_binsize, q = seg_robost_quant)
    print(f'possible bp_pair more than {len(reads_with_SA)}')
    if len(reads_with_SA) > 0:
        sa_bp_pairs =  getbppair.get_bp_pair(reads_with_SA)
    else:
        sa_bp_pairs = pd.DataFrame()

    if bamfile_init.downsample_ratio != 1:
        bamfile = ds_bamfile_pathname

    logg.info(f'process_alignment finished')
    ## return
    return bp_all, clip_bg, bin_count, sa_bp_pairs, bamfile


#------------------------------------------------------------------------------#
def detect_amplified_segment(bin_count, bp_all, blfile = None, genefile = None,
                             gsfile = None, cgfile = None, biasfile = None,
                             binsize = 10000,
                             min_cn = None, std_scale = 8,
                             min_segsize_scale = 5,
                             min_nt = 5, bamfile = None):
    """detect amplified segments and breakpoint list

    Parameters
    ----------
    bin_count : bin counts
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    genefile : gene annotation file in bed12 format
    gsfile : genome size file (not yet in use)
    cgfile : cancer gene annotation file, gene symbol in 1st column
    binsize : genomic bin size to count reads in
    min_cn : minimal copy number to treat bin as amplification
        score based on assuming no whole genome duplication
        (None - use mean + std_scale * std as cutoff)
    std_scale : determin min_cn by mean + std_scale * std
    min_segsize_scale: minimal segment size factor
        min_segsize_scale * binsize
    min_nt : minimal nt required to be matched to bridge two breakend
    """
    logg.info(f'detect_amplified_segment ....')
    # get the bias in binsize and remove bias from count data
    if (biasfile is not None):
        dfm = cndetect.get_bin_bias(bin_count, biasfile, binsize)
        if not dfm.empty:
            bin_norm = cndetect.norm_glm(dfm)
        else:
            bin_norm = bin_count
    else:
        bin_norm = bin_count

    # set blacklist bins to diploid CN
    if (blfile is not None):
        bin_norm = cndetect.bin_filter_by_blacklist(
            bin_norm, blfile, binsize = binsize)

    # smoothing on CN column
    bin_norm = cndetect.smooth_cn(bin_norm)

    # segmentation
    cn_seg = cbs.cbs(bin_norm, binsize = binsize)

    # get amplicon
    cn_amp_raw = cbs.cbs_amplicon(cn_seg, min_cn = min_cn)

    # filter out small spikes (fold of binsize)
    cn_amp_raw = cndetect.segment_filter_by_size(cn_amp_raw, \
        binsize = binsize, fold = min_segsize_scale)

    # cn_seg boundary refine
    cn_amp = cndetect.finemap_cn_segment_boundary(cn_amp_raw, bp_all, \
        binsize = binsize)

    # annotation and output
    if (genefile is not None):
        cn_amp = cndetect.segment_annotation(cn_amp, genefile)
        if (cgfile is not None):
            cn_amp = misc.map_cancergene(cn_amp, cgfile)

    logg.info(f'detect_amplified_segment finished')
    ## return
    return cn_amp, cn_seg, bin_norm


#-----------------------------------------------------------------------------#
# def compute_depth(bin_norm, cn_amp, cn_seg):
#     table_df1 = bin_norm.agg({'Count': ['sum', 'mean'],
#                         'CN': ['sum', 'mean']}).T
#     table_df1.columns = ['bin_norm_sum','bin_norm_mean']
#     table_df2 = cn_amp.agg({'Count': ['sum', 'mean'],
#                             'CN': ['sum', 'mean']}).T
#     table_df2.columns = ['cn_amp_sum','cn_amp_mean']
#     table_df3 = cn_seg.agg({'Count': ['sum', 'mean'],
#                             'CN': ['sum', 'mean']}).T
#     table_df3.columns = ['cn_seg_sum','cn_seg_mean']
#     table_df = pd.concat([table_df1, table_df2, table_df3], axis=1)
#     return table_df

#-----------------------------------------------------------------------------#
def detect_bp_pair(bamfile, sa_bp_pairs, bp_all, cn_amp,
                   blfile = None, bp_min_clip = 10,
                   binsize = 10000, max_n_bp = 500, clean_bp_perc = 0.2,
                   only_keep_clean_bp = True, min_nt = 5,
                   mapq = 1, nmmax = 3,
                   num_chrom = None, num_stratify = None, no_sa = False,
                   multi_process = True, num_cores = 8,
                   pe_window = 150, bin_count=None, bin_norm=None):
    """ from breakpoint list to breakpoint pairs

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    sa_bp_pairs : SA-derived breakpoint-pair candidates from process_alignment
    no_sa : disable SA candidates and SA-guided discordant-edge snapping
        Clip-site pairing and PE/DS evidence remain enabled.
    bp_all : all breakpoints in a dict
        dict{'left' = pd_df, 'right' = pd_df}
    cn_amp : amplified segments in pd_df
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    bp_min_clip : minimal clipped reads of breakpoint to include
    binsize : genomic bin size
    max_n_bp : max # of breakpoints to included in the analysis
        if the sequence depth is too high, the clip points will be too many
        in this case, only the top N breakpoints will be analyzed
        if tie @ N all breakpoints <= counts @ N will be removed
    clean_bp_perc : percentage cutoff of the "CleanBP"
        CleanBP (breakpoint) is defined as the ClipDepth/InDepth >= perc
        - calculate the sequence depth on the breakpoint (InDepth)
        - # of clip reads (ClipDepth)
    only_keep_clean_bp: remove non CleanBP
    min_nt : minimal nt required to be matched to bridge two breakend
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)

    """
    logg.info(f'detect_bp_pair ....')
    # print(cn_amp.iloc[:,:6])
    cn_amp = cn_amp[['Chrom', 'Start', 'End', 'CN', 'ClipLeft', 'ClipRight', 'Gene', 'CancerGene']]
    # bamfile_ = pysam.AlignmentFile(bamfile, 'r')
    # bamfile_init = covbam.prebamfile(bamfile_, downsample=15)
    # bamfile_init.downsample_ratio = downsample_ratio
    print(f'cn_amp:{cn_amp}')

    def estimate_insert_stats(bam, n_sample=1000, seed=0):
        """
        Estimate read length, insert size, and insert-size standard deviation from proper pairs
        """
        rng = np.random.default_rng(seed)
        lengths, inserts = [], []
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped:
                continue
            lengths.append(a.infer_query_length(always=False) or 0)
            if a.is_proper_pair and not a.is_reverse and 0 < a.template_length < 10000:
                inserts.append(a.template_length)
            if len(inserts) >= n_sample:
                break
        read_length = float(np.mean([l for l in lengths if l])) if lengths else 100.0
        insert_size = float(np.mean(inserts)) if inserts else 300.0
        insert_std = float(np.std(inserts)) if inserts else 50.0
        max_insert = insert_size + 3 * insert_std
        return read_length, insert_size, insert_std, max_insert


    regions = []
    for chrom, start, end in cn_amp.iloc[:, :3].itertuples(index=False,name=None):
        start = max(0, int(start) - 10 * 10000)
        end = int(end) + 10 * 10000
        regions.append( (str(chrom), start, end) )

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        read_length, _, _, max_insert = estimate_insert_stats(bam, n_sample = 5000)
        edges = covbam.cluster_discordant_edges_multi(bam, regions = regions, max_insert = max_insert, read_length = read_length,
                                                      pair_support = max(2, bp_min_clip))
        # Large WGS amplicons can produce tens of thousands of PE edges.  Log
        # the count instead of serializing the full list (diagnostic only).
        print(f'discordant edge count: {len(edges)}')

    if not no_sa:
        edges, _ = cndetect.snap_discordant_edges_to_bp_duo(edges, sa_bp_pairs)

    cn_amp, edges = cndetect.refine_boundary_with_discordant( cn_amp, edges, bin_count, binsize = 3*10000,
                                                             min_edge_support = max(2, bp_min_clip), partner_anchor_window =  10000)
    ##----------------------------
    if not no_sa:
        # get bp pairs from bp_cand_stats
        if not sa_bp_pairs.empty:
            print(f'get the number of bppair from bamfile  {sa_bp_pairs.shape}')
            # filter breakpoints in blacklist regions
            if (blfile is not None):
                sa_bp_pairs = bpdetect.bp_filter_by_blacklist(sa_bp_pairs, blfile)
            print(f'get the number of bppair after filter_by_blacklist {sa_bp_pairs.shape}')

            # oxnly include the breakpoints in amplified segments
            sa_bp_pairs = bpdetect.bp_filter_by_amplicon(sa_bp_pairs, cn_amp, binsize)
            print(f'get the number of bppair after filter_by_amplicon {sa_bp_pairs.shape}')

            # filter when the bp pair is very large
            if sa_bp_pairs.shape[0] >= 20000:
                sa_bp_pairs = sa_bp_pairs[sa_bp_pairs["Count"] > 1]
            bp_cand_bam =  bpdetect.refine_bp_cand_bam(sa_bp_pairs, bp_all)  # sa_bp_pairs[sa_bp_pairs["Count"] > 1]
        else:
            bp_cand_bam = pd.DataFrame()
    else:
        sa_bp_pairs, bp_cand_bam = pd.DataFrame(), pd.DataFrame()

    # get the breakpoint candidates >= "n" read counts
    bp_cand = bpdetect.bp_filter_by_n(bp_all, n = bp_min_clip, edges=edges, pe_window = pe_window)
    print(f'bp_cand after filter breakpoints  >= "n" read counts: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')

    # filter breakpoints in blacklist regions
    if (blfile is not None):
        bp_cand = bpdetect.bp_filter_by_blacklist(bp_cand, blfile, isdf = False)
    print(f'bp_cand after filter breakpoints in blacklist regions: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')

    # only include the breakpoints in amplified segments
    bp_cand = bpdetect.bp_filter_by_amplicon(bp_cand, cn_amp, binsize, isdf = False)
    print(f'bp_cand after filter breakpoints in amplified segments: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')

    # in this case, only the top N breakpoints will be analyzed
    bp_cand_df = bpdetect.bp_top_n(bp_cand, topn = max_n_bp, num_chrom = num_chrom, num_stratify = num_stratify)
    print(f'bp_cand_df.shape: {bp_cand_df.shape}')

    # put bp candidates into one table  ----------------------(consuming time )
    bp_cand_all = bpdetect.bp_seq_depth(bp_cand_df, \
                                        bamfile, perc = clean_bp_perc, mapq = mapq, nmmax = nmmax, isdf = False)
    # remove noise
    if (only_keep_clean_bp):
        bp_cand_stats = bp_cand_all[bp_cand_all['CleanBP'] == True]
    else:
        bp_cand_stats = bp_cand_all


    print(f'get bp pairs')
    # get bp pairs
    bp_duo = bpjoint.get_bppair(bamfile, bp_cand_stats, sa_bp_pairs, min_nt = min_nt,
                                multi_process = multi_process, num_cores = num_cores)  # consuming time

    print(f'bp_duo.shape: {bp_duo.shape}')

    # merge bp cand
    try:
        bp_cand_stats.drop(columns = 'Group', inplace=True)
    except:
        print(f'bp_cand_stats no Group column or empty')

    if not bp_cand_bam.empty:
        bp_cand_stats = pd.concat([bp_cand_stats, bp_cand_bam]).drop_duplicates(['Chrom', 'Coord', 'Clip'], keep='first')
    # get bp pairs from paired reads
    bp_cn = bpjoint.bp_cn_boundary(cn_amp)
    bp_duo_pe = bpjoint.get_bppair_peread(bamfile, bp_cn, min_reads = 2)
    
    edges_df = evidence.edges_to_bp_pair_df(edges)
    # Keep split/bp_duo junction metadata when another evidence source reports
    # the same pair. If no bp_duo exists, prefer PE over DS; all support is merged.
    bp_duo = evidence.merge_bp_pair_evidence(bp_duo, bp_duo_pe, edges_df)
    print(f'merged DS/PE/split breakpoint pairs: {bp_duo.shape}')

    # Snap one unsupported DS endpoint to a nearby clip-supported CN boundary.
    if bin_norm is not None:
        bp_duo = evidence.snap_ds_to_cn_boundary(
            bp_duo, bp_cand_stats, cn_amp, bin_norm, bp_all,
            min_cn=4.5, min_reads=2, min_boundary_clip=5,
            max_snap_distance=30000, max_cn_gap=30000,
            min_continuous_cn_frac=1.0)
    else:
        print('DS-to-CN boundary snapping skipped: bin_norm is unavailable')

    logg.info(f'detect_bp_pair finished')
    return bp_duo, bp_cand_stats,cn_amp

# @profile
def construct_amplicon(bp_duo, bp_cand_stats, cn_amp, bin_norm, 
                       genefile = None, sefile = None, cgfile = None, 
                       min_junc_cnt = 10, force_juncread = False, bpp_min_dist = 200,
                       segment_restrict = True, subseg = True, rm_amp_df = False,
                       knn = 3, saveseg = None, savebp_bp_pair = None, EGA = None,
                       min_cn = 5, min_circ_score = 4.0,
                       precomputed_seg = None):
    """construct amplicon structure

    Parameters
    ----------
    bp_duo : final breakpoint pairs used to construct amplicon
    bp_cand_stats: breakpoint candidate stats table
    cn_amp : copy number amplified segments
    genefile : gene annotation file
    cgfile : cancer gene file with first column gene symbol
    min_junc_cnt : minimal supporting read counts of breakpoint pairs
    bpp_min_dist : remove breakpoint pairs with distance <= N
    segment_restrict : True - remove non CleanBP; False - use non CleanBP as alternative segment stop
    min_circ_score : float, default 4.0
        Keep circular amplicons with Score >= this threshold. Set to
        ``-inf`` to disable the final score filter.
    """
    logg.info(f'construct_amplicon ...')
    #--------------------------------------------------------------------------#
    # return empty construct if no bp pairs detected
    if (bp_duo.empty):
        colnames = ['Chrom', 'Start', 'End', 'Strand',
                    'SplitCount', 'CN', 'AmpliconID', 'Gene']
        colnames1 = ['AmpliconID', 'Chrom1','Start', 'Chrom2', 'End', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 
                     'SplitCount_std', 'CN_sum', 'CN_mean', 'CN_std', 
                     'FCleft_sum', 'FCright_sum', 'invCNCV_sum', 'invCNCV_mean', 'invSplitCV',
                     'Gene_num', 'Cancergene_num', 'SE_num', 'FCleft_mean_1', 'FCright_mean_1', 'Score']
        circ_anno = pd.DataFrame(columns = colnames)
        line_anno = pd.DataFrame(columns = colnames)
        circ_score = pd.DataFrame(columns = colnames1)
        return circ_anno, line_anno, bp_duo, cn_amp, circ_score

    #--------------------------------------------------------------------------#
    required_col = [ 'Chrom', 'Start', 'End', 'CN', 'ClipLeft', 'ClipRight', 'Gene', 'CancerGene']
    optional_col = ['DiscLeft', 'DiscRight', 'PartnerLeft', 'PartnerRight', 'SegSource']
    cn_amp = cn_amp[ required_col + [x for x in optional_col if x in cn_amp.columns] ]

    print(f'filter bp duo: bp_duo.shape {bp_duo.shape}')

    # estimate the thred for bp_pair
    # force_juncread: use the caller-supplied min_junc_cnt (--juncread) as-is,
    #   bypassing the auto-estimated / hardcoded threshold below.
    if force_juncread:
        print(f'force_juncread: use min_junc_cnt = {min_junc_cnt} (skip estimaed_thred)')
        aa = bppair.estimaed_thred(bp_duo[bp_duo['Seq'] !='DS'])
        print(f'Estimated threshold: {aa}')
    elif bp_duo.shape[0] >= 500:
        _, min_junc_cnt, _, _ = bppair.estimaed_thred(bp_duo[bp_duo['Seq'] !='DS'])
    else:
        min_junc_cnt = 2
    
    if EGA:
        # keep = bp_duo['Count'] > min_junc_cnt
        # if 'Seq' in bp_duo.columns:
        #     keep &= bp_duo['Seq'].ne('DS')
        # print(f'EGA breakpoint filter only: Count > {min_junc_cnt}, exclude Seq == DS')
        # bp_pair = bp_duo[keep].copy()
        # if 'ClipDepth' in bp_cand_stats.columns:
        #     bp_cand_stats = bp_cand_stats[
        #         bp_cand_stats['ClipDepth'] >= min_junc_cnt].copy()
        #     print(f'EGA breakpoint candidates: ClipDepth >= {min_junc_cnt}')
        # bp_pair = bppair.bp_pair_filter(bp_duo, cn_amp, bp_cand_stats, 0, min_junc_cnt,
        #                                max_srp=1)
        bp_pair = bppair.bp_pair_filter(bp_duo, cn_amp, bp_cand_stats, bpp_min_dist, min_junc_cnt)
    else:
        bp_pair = bppair.bp_pair_filter(bp_duo, cn_amp, bp_cand_stats, bpp_min_dist, min_junc_cnt)
    # bp_pair = bppair.snap_pe_to_segment_boundary(bp_pair, cn_amp, tolerance=20)

    logg.info(f'filter bp_pair by thred: {min_junc_cnt}: bp_pair.shape {bp_pair.shape}')    # if subseg and bpp_max_dist:
    #     print(f'only subseg is {subseg} == True')
    #     ind = (bp_pair['Chrom1']==bp_pair['Chrom2']) & (np.abs(bp_pair['Coord2'] - bp_pair['Coord1']) > bpp_max_dist)
    #     bp_pair = bp_pair[~ind]

    # breakpoint stats table
    logg.info(f'bp_refine: bp_pair.shape {bp_pair.shape}')
    bp_fine, cn_amp_merged = bppair.bp_refine(bp_pair, bp_cand_stats, cn_amp)
    if precomputed_seg is not None:
        seg = precomputed_seg[['Chrom', 'Start', 'End', 'CN']].copy()
        logg.info(f'use precomputed breakpoint segments: {seg.shape}')
    else:
        seg = bppair.get_segment(bp_fine, bp_pair, cn_amp_merged,
                                 restrict=segment_restrict, subseg=subseg,
                                 rm_amp_df=rm_amp_df, knn=knn)
        seg = bppair.add_cn(seg.copy(), cn_amp, bin_norm)
        # if not EGA:
        seg = bppair.filter_seg_by_interior_cn(seg, bin_norm, min_cn=min_cn, min_amp_frac=0.5)
    if saveseg:
        seg.to_csv(saveseg, sep='\t', index=False)
    if savebp_bp_pair:
        bp_pair.to_csv(savebp_bp_pair, sep='\t', index=False)

    #--------------------------------------------------------------------------#
    print(f'get all circles in the ggraph')
    gg, all_circ_path, circ_anno, second_find = ggraph.find_circ(
        bp_pair, bp_fine, seg, genefile, cgfile, sefile)

    # compute fretures of  circs and  take the better circs based on scores
    if not circ_anno.empty:

        circ_score = ggraph.add_stats_circ(circ_anno, bin_norm, binsize = 10000)
    
        # search the circ again on specfic chroms 

        circ_new = pd.DataFrame()
        if second_find:
            result = (circ_score.groupby(['Chrom1','Start','End']).apply(lambda x: x.nlargest(1, 'Score')).reset_index(drop=True))

            select_circ = result['AmpliconID'][(result['Score'] >= result['Score'].quantile(0.85))] #  & (result['Score'] >= 3
            print(f'select_circ: {select_circ}')


            cand_serch_circ = circ_score[circ_score['AmpliconID'].isin(select_circ)][['AmpliconID', 'Chrom1', 'Start', 'End']]
            select_chrom = cand_serch_circ['Chrom1'].tolist()
            
            bp_pair1 = bp_pair[(bp_pair['Chrom1'].isin(select_chrom)) & (bp_pair['Chrom2'].isin(select_chrom))]
            seg1 = seg[seg['Chrom'].isin(select_chrom)]
            bp_fine1 = bp_fine[bp_fine['Chrom'].isin(select_chrom)]
            _, _, circ_anno1, _ = ggraph.find_circ(
                bp_pair1, bp_fine1, seg1, genefile, cgfile, sefile)
            
            # merge the circ of two search
            query_circ_circ = np.unique(circ_anno['AmpliconID'])
            query_circ_list = []
            for i in query_circ_circ:
                query_circ_list_ = []
                for i in circ_anno[ circ_anno['AmpliconID'] == i ].itertuples():
                    query_circ_list_.append(i[1:7])
                query_circ_list.append(query_circ_list_)

            circ_cand = []
            circ_id = np.unique(circ_anno1['AmpliconID'])
            circ_num = len(query_circ_circ)
            for i in circ_id:
                circ_cand_ = []
                circ_df_sub = circ_anno1[ circ_anno1['AmpliconID'] == i ]
                for j in circ_df_sub.itertuples():
                    circ_cand_.append(j[1:7])
                if circ_cand_ not in query_circ_list:
                    circ_reindex = 'circ_' + str(circ_num)
                    circ_df_sub['AmpliconID'] = circ_reindex
                    circ_cand.append( circ_df_sub )
                    circ_num += 1
            if circ_cand:
                circ_new = pd.concat(circ_cand)
                # compute fretures of  circs and  take the better circs based on scores
                circ_new_score = ggraph.add_stats_circ(circ_new, bin_norm, binsize = 10000)
            else:
                circ_new = pd.DataFrame()

        if not circ_new.empty:
            circ_anno = pd.concat([circ_anno, circ_new])
            circ_score = pd.concat([circ_score, circ_new_score])
    else:
        circ_score = pd.DataFrame()

    # Apply one circle-level score filter to both the segment table and its
    # summary table before any downstream output, matrix generation, plotting,
    # or serialization. This prevents a low-score AmpliconID from surviving in
    # one representation but not the other. The boundary is inclusive:
    # Score == min_circ_score is retained.
    if not circ_score.empty:
        score = pd.to_numeric(circ_score['Score'], errors='coerce')
        keep_ids = set(circ_score.loc[score >= min_circ_score, 'AmpliconID'])
        n_before = circ_score['AmpliconID'].nunique()
        circ_score = circ_score[circ_score['AmpliconID'].isin(keep_ids)].copy()
        circ_anno = circ_anno[circ_anno['AmpliconID'].isin(keep_ids)].copy()
        n_after = len(keep_ids)
        print(
            f'filter circular amplicons by Score >= {min_circ_score:g}: '
            f'{n_before} -> {n_after} (removed {n_before - n_after})'
        )

    #--------------------------------------------------------------------------#
    # Cluster scored circles before downstream output, plotting and matrices.
    # Keep the five-value return interface; audit records also survive pickle.
    all_circ_anno, all_circ_score = circ_anno, circ_score
    circ_anno, circ_score, circle_map = ggraph.deduplicate_ask_circles(
        all_circ_anno, all_circ_score, overlap=0.95,
        endpoint_tolerance=None, strict_overlap=True,
    )
    derived_ids = set(circle_map.loc[circle_map['role'] == 'derived', 'AmpliconID'])
    circ_anno.attrs['circle_dedup'] = {
        'overlap': 0.95, 'strict_overlap': True,
        'endpoint_tolerance': None,
        'map_columns': list(circle_map.columns),
        'map_records': circle_map.to_dict('records'),
        'segment_columns': list(all_circ_anno.columns),
        'derived_segments': all_circ_anno.loc[
            all_circ_anno['AmpliconID'].isin(derived_ids)
        ].to_dict('records') if derived_ids else [],
        'score_columns': list(all_circ_score.columns),
        'derived_scores': all_circ_score.loc[
            all_circ_score['AmpliconID'].isin(derived_ids)
        ].to_dict('records') if derived_ids else [],
    }
    print(f'cluster circular amplicons (reciprocal overlap >0.95): '
          f'{len(circle_map)} -> {len(circ_score)} representatives; '
          f'{len(derived_ids)} derivatives')

    # remove redundant
    line_uniq = gg.path_unique(all_circ_path, type = 'linear')
    print(f'get representitive line in the ggraph')
    # choose a representitive one
    line = gg.get_representitive_path(line_uniq, type = 'linear')

    # make interpretable amplicon
    line_df = gg.make_amplicon_df(line, type = 'linear')

    # annotate circular amplicon
    if (genefile is not None):
        line_anno = line_df.assign(
            Gene = grange.GRange.map_gene(line_df, genefile)['Gene'])
        if (cgfile is not None):
            line_anno = misc.map_cancergene(line_anno, cgfile)
    else:
        line_anno = line_df
    if (sefile is not None):
        line_anno = grange.GRange.map_senchancer(line_anno, sefile)

    logg.info(f'construct_amplicon finished')

    return circ_anno, line_anno, bp_pair, seg, circ_score

#------------------------------------------------------------------------------#
def plot_amplicon(circ_anno, line_anno = None, cn_amp = None,
                  genefile = None, sefile = None, bincnt = None, binsize = 10000,
                  ext = 0.3, fig_dir = 'plot', plot_n = 5,
                  fig_width = 15, fontsize = 12):
    """plot amplicon structure

    Parameters
    ----------
    circ_anno : circular amplicon dataframe
    line_anno : linear amplicon dataframe
    cn_amp : amplified segments dataframe
    genefile : gene annotation file
    sefile: super enhancer annotation file *********************  no set
    bincnt : bin read counts per binsize genomic region
    ext : extent roi by 30 percent
    plot_dir : plot folder
    plot_n : only plot n structures for linear or cn
        Note : already plot all structures for circular
    """

    #--------------------------------------------------------------------------#
    # readin bed12 gene annotation file
    if (genefile is not None):
        genebed = trackplot.read_bedf(genefile)
        genebed = [g for g in genebed if (g.name != '')]
    
    # readin bed4 super encancer annotation file
    if (sefile is not None):
        sebed = trackplot.read_bedf(sefile)
        sebed = [g for g in sebed if (g.name != '')]

    #--------------------------------------------------------------------------#
    # plot amplified segments
    if (cn_amp is not None and not cn_amp.empty):
        # sort cn_amp by copy number gain
        cn_amp_list = trackplot.get_rois_with_score(cn_amp, binsize * 10)
        cn_amp_list = sorted(cn_amp_list, key=lambda x: x[3], reverse = True)

        # plot
        for tag, roi in enumerate(cn_amp_list[:plot_n]):
            # get sub dataframe
            df_sub = trackplot.get_sub_df(cn_amp, roi)

            # output pdf file name
            fig_out = os.path.join(fig_dir, 'ampseg_' + str(tag) + '.pdf')

            # prepare plot
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                df_sub, genebed, bincnt, binsize, ext, 'cn')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot circular amplicon
    if (not circ_anno.empty):
        print(circ_anno)
        for tag, circ_sub in circ_anno.groupby(['AmpliconID']):

            # output pdf file name
            try:
                fig_out = os.path.join(fig_dir, 'circular_' + tag + '.pdf')
            except:
                fig_out = os.path.join(fig_dir, 'circular_' + tag[0] + '.pdf')
            # prepare plot
            circ_sub_copy = circ_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                circ_sub_copy, genebed, bincnt, binsize, ext, 'circular')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot linear amplicon
    if (line_anno is not None and not line_anno.empty):
        id2plot = list(
            line_anno.groupby('AmpliconID').CN.max().\
                sort_values(ascending = False)[0:plot_n].index)
        for tag in id2plot:
            # output pdf file name
            try:
                fig_out = os.path.join(fig_dir, 'Linear_' + tag + '.pdf')
            except:
                fig_out = os.path.join(fig_dir, 'Linear_' + tag[0] + '.pdf')
            # get the sub dataframe
            line_sub = line_anno[line_anno['AmpliconID'] == tag]

            # prepare plot
            line_sub_copy = line_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                line_sub_copy, genebed, bincnt, binsize, ext, 'linear')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)


if __name__ == '__main__':
    datapath = '/Users/weina/ask-main/data'
    resultpath = '/Users/weina/project/ecDNA/ask_run/test_result/GBM39_chr8_all_ask'
    resfigpath = '/Users/weina/project/ecDNA/ask_run/test_result/GBM39_chr8_all_ask'
    from pathlib import Path
    import os
    import time
    from pathlib import Path
    import pickle
    #------------------------------------------------------------------------#
    #@ create base dir of output files
    Path(os.path.dirname(resultpath)).mkdir(parents=True, exist_ok=True)

    #@ create base dir of plot files
    Path(resfigpath).mkdir(parents=True, exist_ok=True)

    #@ set output file names
    output_stats     = resultpath + '_ask_stats.tsv'
    output_cnamp     = resultpath + '_ask_amplified_segment.tsv'
    output_cnseg     = resultpath + '_ask_cn_segmentation.tsv'

    output_bpcand    = resultpath + '_ask_breakpoint.tsv'
    output_bpduo     = resultpath + '_ask_breakpoint_pair_raw.tsv'
    output_bppair    = resultpath + '_ask_breakpoint_pair.tsv'
    
    output_circ      = resultpath + '_ask_amplicon_circular.tsv'
    output_circ_stat = resultpath + '_ask_amplicon_circular_stat.tsv'
    output_line      = resultpath + '_ask_amplicon_linear.tsv'
    output_seg      = resultpath + '_ask_amplicon_seg.tsv'
    output_clip      = resultpath + '_ask_clip_count.bedgraph'
    output_bincnt    = resultpath + '_ask_bin_count.tsv'
    output_binnorm   = resultpath + '_ask_bin_count_norm.tsv'
    output_align     = resultpath + '_ask_breakpoint_pair_alignment.tsv'

    output_pdat_1    = resultpath + '_ask_step1.pdat'
    output_pdat_2    = resultpath + '_ask_step2.pdat'
    output_pdat_3    = resultpath + '_ask_step3.pdat'
    output_pdat_4    = resultpath + '_ask_step4.pdat'

    
    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#

    # blacklist file
    blfile = os.path.join(datapath,  'hg38_blacklist.bed')

    # gene annotation file
    genefile = os.path.join(
        datapath, 'hg38_refgene.bed12.gz')

    # genome size file (currently not in use)
    gsfile = os.path.join(
        datapath, 'hg38.genome')

    # cancer gene file
    cgfile = os.path.join(
        datapath, 'Census_all_20200624_14_22_39.tsv')

    bamfile = os.path.join(
        datapath, 'test_sort.bam')

    # bias file (GC + mappability)
    biasfile = os.path.join(
         datapath,
        'hg38_bias.bed.gz')

    binsize = 10000
    mapq = 20
    nmmax = 1
    segmode = 'standard'
    sub_binsize = 1000
    segquant = 0.5
    mincn = 5
    bpcount = 5
    juncread = 5
    max_n_breakpoint = 500
    cleanbp_cutoff = 0.2
    sub_binsize = 1000
    only_keep_cleanbp = True

    startTime = time.time()
    #--------------------------------------------------------------------------#
    # step 1 : process bam file to generate bin counts and clip reads
    #--------------------------------------------------------------------------#
    # run module
    bp_all, clip_bg, bin_count, reads_with_SA = process_alignment(
        bamfile, gsfile = gsfile, binsize = binsize,
        mapq = mapq, nmmax = nmmax,
        mode = segmode, sub_binsize = sub_binsize,
        seg_robost_quant = segquant, SA_with_nm=True)

       
    # output results
    bin_count.to_csv(output_bincnt, sep='\t', index=False)
    clip_bg.to_csv(output_clip, sep='\t', index=False, header=False)

    # save pdat
    pdat = [bp_all, clip_bg, bin_count]
    with open(output_pdat_1, "wb") as f:
        pickle.dump(pdat, f)

    # output log
    print('alignment process - done \
        - {0:0.1f} seconds\n'.format(time.time() - startTime))

    with open('GBM39_chr78_ask_0927/GBM39_chr78_ask_ask_step1.pdat', 'rb') as f:
        dat = pickle.load(f)
    bp_all, clip_bg, bin_count = dat[0],dat[1],dat[2]
    #--------------------------------------------------------------------------#
    # step 2 : identify amplified segments
    #--------------------------------------------------------------------------#
    # run module
    cn_amp, cn_seg, bin_norm = \
        detect_amplified_segment(
            bin_count, bp_all,
            blfile, genefile, gsfile, cgfile,
            biasfile = biasfile,
            binsize = binsize, min_cn = mincn)

    # output results
    cn_amp.to_csv(output_cnamp, sep='\t', index=False)
    cn_seg.to_csv(output_cnseg, sep='\t', index=False)
    bin_norm.to_csv(output_binnorm, sep='\t', index=False)

    # save pdat
    pdat = [cn_amp, cn_seg, bin_norm]
    with open(output_pdat_2, "wb") as f:
        pickle.dump(pdat, f)

    # output log
    print('detect amplified segments - done \
        - {0:0.1f} seconds\n'.format(time.time() - startTime))

    with open('GBM39_chr78_ask_0927/GBM39_chr78_ask_ask_step2.pdat', 'rb') as f:
        dat = pickle.load(f)
    cn_amp, cn_seg, bin_norm = dat[0],dat[1],dat[2]

    #--------------------------------------------------------------------------#
    # step 3 : identify breakpoint pairs
    #--------------------------------------------------------------------------#
    # run module
    bp_duo, bp_cand_stats = detect_bp_pair(
        bamfile,reads_with_SA, cn_amp, blfile,
        binsize = binsize, bp_min_clip = bpcount,
        mapq = mapq, nmmax = nmmax,
        max_n_bp = max_n_breakpoint,
        clean_bp_perc = cleanbp_cutoff,
        only_keep_clean_bp = only_keep_cleanbp, 
        no_sa = True)

    # output results
    bp_duo.to_csv(output_bpduo, sep='\t', index=False)

    # save pdat
    pdat = [bp_duo, bp_cand_stats]
    with open('both_chr78_ask_step3.pdat', "wb") as f:
        pickle.dump(pdat, f)

    # output log
    print('detect breakpoint pairs - done \
        - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # step 4 : construct amplicons
    #--------------------------------------------------------------------------#
    # run module
    circ_anno, line_anno, bp_pair, seg = construct_amplicon(
        bp_duo, bp_cand_stats, cn_amp, bin_norm, genefile, cgfile,
        segment_restrict = False,
        min_junc_cnt = juncread)
    

    # output results
    circ_anno.to_csv(output_circ, sep='\t', index=False)
    line_anno.to_csv(output_line, sep='\t', index=False)
    bp_pair.to_csv(output_bppair, sep='\t', index=False)
    seg.to_csv(output_seg, sep='\t', index=False)


    # save pdat
    pdat = [circ_anno, line_anno, bp_pair, seg]
    with open(output_pdat_4, "wb") as f:
        pickle.dump(pdat, f)

    # output log
    print('construct amplicons - done \
        - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # step 5 : plot amplicons
    #--------------------------------------------------------------------------#
    # run module
    outfig = resfigpath
    try:
        plot_amplicon(circ_anno, line_anno, cn_amp, genefile,
                          bin_norm, binsize = binsize,
                          fig_dir = outfig, plot_n = 5,
                          fig_width = 15, fontsize = 12, ext = 0.3)
        # output log
        print('plot amplicons - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))
    except:
        print('plot amplicons - error \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # output stats
    #--------------------------------------------------------------------------#
    # output statistics table
    with open(output_stats, 'w') as f:
        f.write('# of amplified segments: ' + str(len(cn_amp)) + '\n')
        f.write('# of candidate breakpoints: ' + str(len(bp_cand_stats)) + '\n')
        f.write('# of breakpoint pairs: ' + str(len(bp_duo)) + '\n')
        f.write('# of circular amplicons: '
            + str(circ_anno['AmpliconID'].nunique()) + '\n')
        f.write('# of linear amplicons: '
            + str(line_anno['AmpliconID'].nunique()) + '\n')

    print('All - done - {} seconds\n'.format(time.time() - startTime))
