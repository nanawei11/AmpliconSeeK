################################################################################
# infer segments from breakpoint pairs
################################################################################
#------------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import pysam
import math
import scipy.stats as stats
import statsmodels.api as sm
from itertools import groupby
from itertools import combinations
from operator import itemgetter
from collections import defaultdict
from statsmodels.formula.api import glm
import matplotlib.pyplot as plt
from collections import defaultdict

import cbs
import misc
import covbam
from grange import GRange
from bpjoint import join_breakpoint


#------------------------------------------------------------------------------#
# get bias in binsize
#------------------------------------------------------------------------------#
def get_bin_bias(bin_count, biasfile, binsize = 10000):
    """calculate bias factor in binsize
    """
    df = pd.read_csv(biasfile, sep = '\t')
    df.columns = ['Chrom', 'Start', 'End', 'GCcontent', 'Mappability']
    df.Start = (np.floor(df.Start/binsize)*binsize).astype('int')
    df_agg = df.groupby(['Chrom', 'Start'])[['GCcontent', 'Mappability']]\
        .mean().reset_index().rename(columns = {"Start":"Coord"})
    dfm = pd.merge(bin_count, df_agg, on = ['Chrom', 'Coord'])
    return dfm


#------------------------------------------------------------------------------#
# glm remove GC and mappability biasa
#------------------------------------------------------------------------------#
def norm_glm(dfm, cutoff = [0, 0], method = 'nb', fold_bound = 5, plot = False):
    """fit glm to remove biases from bin counts

    Parameters
    ----------
    cutoff : select data points used in the model, default use all
    method : 'nb' - negative binomial
             'poisson' - poisson
    fold_bound : bound of bias factor to consider
        to avoid getting super small or large normalized counts
    """
    keep = (dfm.GCcontent >= cutoff[0]) & (dfm.Mappability >= cutoff[1])
    if (method == 'nb'):
        fit = glm('Count ~ GCcontent + Mappability', data = dfm,
            subset = keep, family = sm.families.NegativeBinomial()).fit()
    else:
        fit = glm('Count ~ GCcontent + Mappability', data = dfm,
            subset = keep, family = sm.families.Poisson()).fit()

    coef = fit.params
    fitted_y = np.exp(coef[0] + coef[1] * dfm.GCcontent +
                      coef[2] * dfm.Mappability)
    fac = fitted_y / np.median(fitted_y)

    fac[fac > fold_bound] = fold_bound
    fac[fac < 1/fold_bound] = 1/fold_bound

    norm_y = round(dfm.Count / fac, 4)

    bin_norm = dfm.iloc[:,0:4].copy()
    bin_norm['Count'] = norm_y + 1
    
    norm_fac = get_mode(bin_norm['Count'])
    if norm_fac == 0:
        norm_fac = np.mean(bin_norm['Count'])
    bin_norm['CN'] = 2 * bin_norm['Count'] / norm_fac

    if (plot):
        plt.plot(dfm.GCcontent, dfm.CN, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.GCcontent, norm_y, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.Mappability, dfm.CN, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.Mappability, norm_y, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

    return bin_norm


#------------------------------------------------------------------------------#
# 2D bias plot
#------------------------------------------------------------------------------#
def plot_bias_2D(x, y, z, npts = 100, cmap = 'gist_rainbow',
                 xlab = 'GCcontent', ylab = 'Mappability'):
    """ 2D bias plot ('viridis')
    """

    # bin data and take average
    zi = stats.binned_statistic_2d(y, x, z, 'mean', bins=npts).statistic
    print(np.nansum(abs(zi-np.nanmedian(zi))))

    # make heatmap
    plt.imshow(zi, cmap=cmap, origin = 'lower', extent=(0,1,0,1))
    plt.colorbar()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.show()


#------------------------------------------------------------------------------#
# apply triangular smoothing on the normalized bin counts
#------------------------------------------------------------------------------#
def smooth_cn(bin_norm, k = 10):
    """smooth normalized bin counts
    """
    df = bin_norm.copy()
    # calculate log2 ratio
    df['Log2Ratio'] = np.log2(np.array(df.CN / 2))

    df_smooth = []
    for chr in misc.unique(df['Chrom']):
        dfsub = df[df['Chrom'] == chr]
        x = np.array(dfsub['Log2Ratio'])
        x_smooth = trismooth_1d(x, k = k)
        dfsub = dfsub.assign(Log2Ratio = x_smooth)
        df_smooth.append(dfsub)
    return pd.concat(df_smooth)


#------------------------------------------------------------------------------#
def trismooth_1d(x, k = 10, outlier_sd = 5,
                 edge_repeat = False):
    """triangular smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.asarray(x, dtype=float).copy()
    n = len(x)
    if n == 0 or k <= 0:
        return x

    # Compute neighbour statistics with 1-D convolutions.  The old code
    # materialised an N x (2*k+1) matrix (N x 21 at the default); this keeps
    # only O(N) vectors while preserving the same window and excluding the
    # centre bin from the mean and standard deviation.
    finite = np.isfinite(x)
    values = np.where(finite, x, 0.0)
    neighbour_kernel = np.ones(2 * k + 1, dtype=float)
    neighbour_kernel[k] = 0.0

    if edge_repeat:
        value_pad = np.pad(values, (k, k), mode='edge')
        valid_pad = np.pad(finite.astype(float), (k, k), mode='edge')
    else:
        value_pad = np.pad(values, (k, k), mode='constant')
        valid_pad = np.pad(finite.astype(float), (k, k), mode='constant')

    neighbour_n = np.convolve(valid_pad, neighbour_kernel, mode='valid')
    neighbour_sum = np.convolve(value_pad, neighbour_kernel, mode='valid')
    neighbour_sq_sum = np.convolve(
        value_pad * value_pad, neighbour_kernel, mode='valid')

    ys_mean = np.full(n, np.nan, dtype=float)
    np.divide(neighbour_sum, neighbour_n, out=ys_mean, where=neighbour_n > 0)

    # Sample variance (ddof=1), matching np.nanstd(..., ddof=1).
    ys_var = np.full(n, np.nan, dtype=float)
    valid_var = neighbour_n > 1
    variance_numerator = (
        neighbour_sq_sum - neighbour_sum * neighbour_sum /
        np.maximum(neighbour_n, 1.0)
    )
    np.divide(
        np.maximum(variance_numerator, 0.0), neighbour_n - 1.0,
        out=ys_var, where=valid_var)
    ys_sd = np.sqrt(ys_var)
    ys_upper = ys_mean + outlier_sd * ys_sd
    ys_lower = ys_mean - outlier_sd * ys_sd

    # calculate on outlier
    outlier = finite & ((x > ys_upper) | (x < ys_lower))

    # The centre weight is zero because the outlier itself is excluded before
    # replacement.  Padding followed by valid convolution produces exactly N
    # results even for chromosomes shorter than the window.
    triangular_kernel = np.r_[
        np.arange(1, k + 1, dtype=float),
        0.0,
        np.arange(k, 0, -1, dtype=float),
    ]
    weighted_sum = np.convolve(value_pad, triangular_kernel, mode='valid')
    weight_sum = np.convolve(valid_pad, triangular_kernel, mode='valid')
    avg = np.full(n, np.nan, dtype=float)
    np.divide(weighted_sum, weight_sum, out=avg, where=weight_sum > 0)

    # assign back
    smooth_x = x.copy()
    smooth_x[outlier] = avg[outlier]

    return smooth_x


#------------------------------------------------------------------------------#
def trismooth_1d_plus(x, k = 10, outlier_sd = 5,
                      trim_quantile = 0.1, edge_repeat = False):
    """triangular smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    trim_quantile : remove trim_quantile extreme data points
                    to calcuate mean and sd
                    for both upper and lower quantile
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.array(x)
    y = np.zeros ((len (x), 2*k + 1), dtype = x.dtype)
    y[:, k] = x
    for i in range(k):
        j = k - i
        y[j:,i] = x[:-j]
        y[:-j,-(i+1)] = x[j:]
        if (edge_repeat):
            y[:j,i] = x[0]
            y[-j:,-(i+1)] = x[-1]
        else:
            y[:j,i] = np.nan
            y[-j:,-(i+1)] = np.nan

    # weights
    tt = list(range(1, k+1))
    w = np.array(tt +[k+1] + tt[::-1])

    # smooth
    smooth_x = []
    for i in range(len(x)):
        z = y[i,:]
        zm = z[k]
        z[k] = np.nan

        with np.errstate(invalid='ignore'):
            tf = (z < np.nanquantile(z, trim_quantile)) |\
                 (z > np.nanquantile(z, 1 - trim_quantile))
        z[tf] = np.nan

        z_mean = np.nanmean(z)
        z_sd = np.nanstd(z, ddof=1)
        z_upper = z_mean + outlier_sd * z_sd
        z_lower = z_mean - outlier_sd * z_sd

        if (zm > z_upper or zm < z_lower):
            idx = np.where(~np.isnan(z))[0]
            smooth_x.append(np.average(z[idx], weights = w[idx]))
        else:
            smooth_x.append(zm)

    return smooth_x


#------------------------------------------------------------------------------#
def get_mode(x):
    """get mode from density
    """
    x = np.array(x)
    x = x[np.isfinite(x)]
    density = stats.kde.gaussian_kde(x)
    v = np.linspace(0, np.quantile(x, 0.99), 500)
    density = density(v)
    return v[density.argmax()]


#------------------------------------------------------------------------------#
# output function
#------------------------------------------------------------------------------#
def get_cn_segment(bin_count_clean, binsize = 10000, nperm = 1000, p = 0.01,
        method = 'cutoff', min_cn = None, std_scale = 8):
    """ get copy number segmentations bin count data
    output a dataframe of amplicons

    Parameters
    ----------
    df : counts in genomic bins
    method : method used to perform segmentation
        'cutoff' : use hard cutoff of min_cn to choose amplified regions
        'cbs' : use cbs method to choose amplified regions

    std_scale : mean + std_scale * std as min_cn
    binsize : target genomic bin size to calculate read counts

    output
    ----------
    segmented regions in dataframe
    """
    df = bin_count_clean.copy()

    # set min_cn by mean + n * std
    if (min_cn is None):
        xx = np.array(df['CN'])
        min_cn = np.mean(xx) + std_scale * np.std(xx, ddof = 1)
        print(min_cn)

    # calculate log2 ratio if use cbs
    if (method == 'cbs'):
        df['L2R'] = np.log2(np.array(df.CN / 2))

    # call copy number segmentations in each chromosome
    cnsegs = []
    for chr in set(df['Chrom']):
        # print(chr)
        dfsub = df[df['Chrom'] == chr]
        # seg = []
        if (method == 'cutoff'):
            x = np.array(dfsub['CN'])
            seg = simple_amplicon_segment(x, min_cn)
            seg = [list(dfsub.iloc[i[0],0:2]) +
                   list(dfsub.iloc[i[1]-1,1:2]+binsize) +
                   [i[2]] for i in seg]
        else:
            x = np.array(dfsub['L2R'])
            seg = cbs.cbs_segment(x, nperm=nperm, p=p)
            seg = cbs.cbs_recheck(x, seg, nperm=nperm, p=p)
            # cbs.cbs_plot_segment(x, seg) # plot function
            seg = [list(dfsub.iloc[i[0],0:2]) +
                   list(dfsub.iloc[i[1]-1,1:2]+binsize) +
                   [2**i[2]*2] for i in seg]
        cnsegs.append(pd.DataFrame(seg,
            columns=['Chrom', 'Start', 'End', 'CN']))
        seg_df = pd.concat(cnsegs)

    return seg_df[seg_df.CN >= min_cn]



#------------------------------------------------------------------------------#
# simple segmentation
#------------------------------------------------------------------------------#
def simple_seg(bin_count, binsize = 10000, min_cn = 5, std_scale = 8):
    """ get copy number segmentations bin count data
    output a dataframe of amplicons
    use hard cutoff of min_cn to choose amplified regions

    Parameters
    ----------
    bin_count : counts in genomic bins
    std_scale : mean + std_scale * std as min_cn
    binsize : target genomic bin size to calculate read counts
    min_cn : minimal copy number to consider

    output
    ----------
    segmented regions in dataframe
    """
    df = bin_count.copy()

    # set min_cn by mean + n * std
    if (min_cn is None):
        xx = np.array(df['CN'])
        min_cn = np.mean(xx) + std_scale * np.std(xx, ddof = 1)
        print(min_cn)

    # call copy number segmentations in each chromosome
    cnsegs = []
    for chr in set(df['Chrom']):
        # print(chr)
        dfsub = df[df['Chrom'] == chr]

        # segmentation
        x = np.array(dfsub['CN'])
        seg = simple_amplicon_segment(x, min_cn)
        seg = [list(dfsub.iloc[i[0],0:2]) +
               list(dfsub.iloc[i[1]-1,1:2]+binsize) +
               [i[2]] for i in seg]

        cnsegs.append(pd.DataFrame(seg,
            columns=['Chrom', 'Start', 'End', 'CN']))
        seg_df = pd.concat(cnsegs)

    return seg_df[seg_df.CN >= min_cn]


#------------------------------------------------------------------------------#
# simple calling of amplicons
#------------------------------------------------------------------------------#
def simple_amplicon_segment(x, min_cn = 5, min_dist_to_merge = 2):
    """
    Simple way to identify amplified segments
    score of genomic bins higher than a certain score
    let's say bins with more than 5 copies

    min_dist_to_merge: merge closeby intervals
    """

    xx = np.where(x > min_cn)[0].tolist()
    segments =[]

    # get the continuous numbers in a list with a tolerance
    grp = groupby(zip(xx, xx[1:] + xx[:0]), \
        lambda x: x[1] - x[0] <= min_dist_to_merge)
    for k,g in grp:
        if (k == True):
            gg = [i for i in g]
            segments.append([gg[0][0], gg[-1][1]+1])

    # get the segmentation mean table
    seg = []
    for j in segments:
        seg_start = j[0]
        seg_end = j[1]
        seg_mean = np.mean(x[seg_start:seg_end])
        seg.append([seg_start, seg_end, seg_mean])

    return seg


#------------------------------------------------------------------------------#
def segment_filter_by_size(cn_amp, binsize = 10000, fold = 5):
    """
    filter copy number segments by size (remove small spikes)

    """
    # old version
    # return cn_amp[cn_amp['End'] - cn_amp['Start'] >= fold * binsize]

    cn_amp_merged = misc.merge_bed(cn_amp, gap = 100000)
    cn_amp_drop = pd.DataFrame([
        row for row in cn_amp_merged if (row[2] - row[1] < fold * binsize
        )], columns = cn_amp.columns[0:3])
    df = pd.merge(cn_amp, cn_amp_drop, indicator = True, how = 'left'
        ).query('_merge == "left_only"').drop('_merge', axis = 1)
    return df


#------------------------------------------------------------------------------#
def get_discordant_protected_bins(bamfile, blacklistfile, bin_norm,
                                   binsize=10000, pair_support=5,
                                   max_insert=900, read_length=100,
                                   protect_extend=100000, min_cn_signal=3):
    """
    Find blacklist bins supported by both CN and discordant-pair evidence.

    Blacklists suppress artifacts from repetitive or low-mappability regions,
    but forcing every overlapping bin to diploid CN can remove a true ecDNA
    segment when a junction connects the region to another locus. A bin is
    protected only when independent discordant pairs support the junction and
    its corrected CN before blacklist filtering is at least min_cn_signal.

    Search only inside blacklist intervals and return protected (Chrom, Coord)
    bins that should not be reset to diploid CN.
    """
    gr = GRange(blacklistfile, 'bedfile')
    bl_regions = [(_gr[0], _gr[1].start, _gr[1].stop) for _gr in gr.gr]
    if not bl_regions:
        return set()

    # process_alignment returns a BAM path, not an open pysam object.
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        edges = covbam.cluster_discordant_edges_multi(
            bam, bl_regions, max_insert=max_insert,
            read_length=read_length, pair_support=pair_support)

    protect = set()
    for endA, endB, n, *_ in edges:
        for chrom, pos, _strand in (endA, endB):
            sub = bin_norm[(bin_norm['Chrom'] == chrom) &
                           (bin_norm['Coord'] >= pos - protect_extend) &
                           (bin_norm['Coord'] <= pos + protect_extend) &
                           (bin_norm['CN'] >= min_cn_signal)]
            for coord in sub['Coord']:
                protect.add((chrom, coord))
    return protect


#------------------------------------------------------------------------------#
def bin_filter_by_blacklist(bin_count, blacklistfile, \
    binsize = 10000, ext = 0, set_value = 2):
    """
    set the bin CN to 2 (diploid), if overlap with blacklist

    ext - extend both end on ext * binsize
    """
    bin_count_ = bin_count.copy()

    # convert bed file to gr
    gr = GRange(blacklistfile, 'bedfile')

    # init blacklist bins
    bl_bin = defaultdict(lambda: False)

    # save blacklist bins to dict
    for _gr in gr.gr:
        _start = math.floor(_gr[1].start/binsize) - ext
        _end = math.floor(_gr[1].stop/binsize) + 1 + ext
        for i in range(_start, _end):
            bl_bin[(_gr[0], i*binsize)] = True


    # get the bool vector of blacklist bins
    tf = [bl_bin[row[1:3]]
          for row in bin_count_.itertuples()]

    # set blacklist bins to 0
    # bin_count_.loc[tf, 'Count'] = set_value
    # bin_count_['CN'] = 2*bin_count_['Count']/np.mean(bin_count_['Count'])

    # set blacklist bins CN to 2, left count unchanged
    bin_count_.loc[tf, 'CN'] = set_value

    return bin_count_

#------------------------------------------------------------------------------#
# identify cn segment boundaries by using clip reads
#------------------------------------------------------------------------------#
def finemap_cn_segment_boundary(cn_amp, bp_all, binsize = 10000):
    """
    search break points for single base resolution
    boundaries of amplified segments
    cn_amp: copy number high amplicons
    bp_all: bp_all
    """
    cn_amp = cn_amp[['Chrom', 'Start', 'End', 'Count', 'CN']]  # 5 col
    seg_bdry = []
    for idx, seg in cn_amp.iterrows():
        bdry = list(seg) + [0, 0]

        # left boundary
        left = bp_all['left']
        left = left[left['Chrom'] == seg[0]]
        left = left.assign(dist = abs(seg[1]-left['Coord']))
        cand = left[left['dist'] <= binsize]
        if (not cand.empty):
            cand = cand.sort_values(['dist'])
            cand = cand.sort_values('Count', ascending=False)
            bdry[1] = list(cand['Coord'])[0]
            bdry[5] = list(cand['Count'])[0]  # adj +1

        # right boundary
        right = bp_all['right']
        right = right[right['Chrom'] == seg[0]]
        right = right.assign(dist = abs(seg[2]-right['Coord']))
        cand = right[right['dist'] <= binsize]
        if (not cand.empty):
            cand = cand.sort_values(['dist'])
            cand = cand.sort_values('Count', ascending=False)
            bdry[2] = list(cand['Coord'])[0]
            bdry[6] = list(cand['Count'])[0]  # adj +1

        # output
        seg_bdry.append(bdry)

    colnames = ['Chrom', 'Start', 'End', 'Count', 'CN', 'ClipLeft', 'ClipRight']
    return pd.DataFrame(seg_bdry, columns=colnames)


#------------------------------------------------------------------------------#
def get_bppair_from_cnamp(cn_amp, bamfile, min_nt = 4):
    """
    get the bppairs from the refined cnamp
    as addition to the bppairs identifed by bp_cand
    """
    op = []
    for row in cn_amp.itertuples():
        op.append([row[1], row[2], 'L'])
        op.append([row[1], row[3], 'R'])
    op

    tt = [pair[0] + pair[1] + list(
            join_breakpoint(bamfile, *pair, min_nt = min_nt)) \
        for pair in combinations(op, 2)]

    t2 = [row[0:6] + [row[6][1] + row[7][1]] \
          for row in tt \
              if row[6][1] > 0 and row[7][1] > 0 \
              and row[6][0] == row[7][0]]

    colnames = ["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", 'Count']
    bp_pair_extra = pd.DataFrame(t2, columns = colnames)

    return bp_pair_extra


#------------------------------------------------------------------------------#
def segment_annotation(cn_amp, refgenef):
    """
    annotate the copy number segments by gene
    """
    return GRange.map_gene(cn_amp, refgenef)


#------------------------------------------------------------------------------#
# filter copy number segments by blacklist intervals (currently not in use)
#------------------------------------------------------------------------------#
def segment_filter_by_blacklist(cn_amp, blacklistfile, f_sort=True):
    """
    filter copy number segments by blacklist intervals

    """
    def filter_blacklist(df, blacklistfile):
        gr1 = GRange(df, 'dataframe_hasend')
        gr2 = GRange(blacklistfile, 'bedfile')
        # extend 100 bp on both end of the blacklist
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 100, invert = True)

        return pd.DataFrame(\
            [[row[0]] + [row[1].start, row[1].stop] + list(row[2]) for row in gr.gr]\
            , columns=df.columns)

    df = filter_blacklist(cn_amp, blacklistfile)

    if (f_sort):
        df = df.sort_values(['Chrom', 'Start', 'End'])

    return df

#------------------------------------------------------------------------------#
# change
#------------------------------------------------------------------------------#


def refine_boundary_with_discordant(
        cn_amp,
        edges,
        bin_count,
        binsize=30000,
        min_edge_support=2,
        partner_anchor_window=10000,
        orphan_window=10000,
        edge_merge_window=150,
        min_cn=4.5,
        min_segment_size=100,
        max_segment_size=100000000,
        cn_binsize=10000,
        weak_clip_max=1,
        weak_clip_boundary_window=10000,
        weak_clip_partner_window=150,
        weak_clip_min_edge_support=5):
    """
    Perform three PE-aware operations:

    1. Refine boundaries whose clip count is exactly zero.
    2. Refine positive but weak boundaries when both ClipLeft and ClipRight
       are <= weak_clip_max.
    3. Rescue missing amplified segments using paired orphan L/R endpoints.

    Weak-clip refinement requires:

    - both segment ends have clip <= weak_clip_max;
    - the current side has a positive weak clip;
    - PE endpoint direction matches the segment boundary;
    - PE endpoint is within weak_clip_boundary_window;
    - PE support reaches weak_clip_min_edge_support;
    - PE endpoint extends the segment outward;
    - the PE partner anchors another amplified segment boundary;
    - the newly included interval has median raw CN >= min_cn;
    - the adjusted segment length is reasonable.

    The weak-clip pass and PE-segment rescue share the same raw CN index.
    """
    out = cn_amp.copy().reset_index(drop=True)
    edges = list(edges) if edges is not None else []

    # ------------------------------------------------------------------
    # Normalize segment columns
    # ------------------------------------------------------------------
    for col in ['DiscLeft', 'DiscRight']:
        if col not in out.columns:
            out[col] = 0

    for col in ['PartnerLeft', 'PartnerRight']:
        if col not in out.columns:
            out[col] = None

    if 'SegSource' not in out.columns:
        out['SegSource'] = 'clip'
    else:
        out['SegSource'] = out['SegSource'].fillna('clip')

    # Make repeated calls idempotent. PE-derived rows are reconstructed.
    out = out[out['SegSource'] == 'clip'].copy().reset_index(drop=True)

    # ------------------------------------------------------------------
    # Normalize PE endpoints and index them by chromosome and direction
    # ------------------------------------------------------------------
    edge_end = []
    edge_index = defaultdict(list)

    for edge_id, edge in enumerate(edges):
        if len(edge) < 3:
            continue

        endpoint_a, endpoint_b = edge[0], edge[1]
        support = int(edge[2])

        if support < min_edge_support:
            continue

        for end_id, (this_end, partner) in enumerate(((endpoint_a, endpoint_b), (endpoint_b, endpoint_a))):

            chrom, coord, strand = this_end
            chrom = str(chrom)
            coord = int(coord)
            side = 'R' if strand == '+' else 'L'

            record = (chrom, coord, side, support, partner, edge_id, end_id)

            edge_end.append(record)
            edge_index[(chrom, side)].append(record)

    # ------------------------------------------------------------------
    # Build the raw, unsmoothed CN index once
    # ------------------------------------------------------------------
    cn_index = {}

    if bin_count is not None and not bin_count.empty:
        for chrom, df in bin_count.groupby('Chrom', sort=False):
            coord = pd.to_numeric(df['Coord'], errors='coerce').to_numpy(dtype=float)

            cn_value = pd.to_numeric(df['CN'], errors='coerce').to_numpy(dtype=float)

            if 'Count' in df.columns:
                count_value = pd.to_numeric(df['Count'], errors='coerce').to_numpy(dtype=float)
            else:
                count_value = np.full(df.shape[0], np.nan)

            valid_coord = np.isfinite(coord)
            coord = coord[valid_coord].astype(np.int64)
            cn_value = cn_value[valid_coord]
            count_value = count_value[valid_coord]

            order = np.argsort(coord)

            cn_index[str(chrom)] = (coord[order], cn_value[order], count_value[order])

    def interval_stats(chrom, start, end):
        """
        Return median CN and mean count for raw bins overlapping [start, end).
        """
        if int(end) <= int(start):
            return np.nan, np.nan

        data = cn_index.get(str(chrom))
        if data is None:
            return np.nan, np.nan

        pos, cn_value, count_value = data

        # A bin [pos, pos + cn_binsize) overlaps [start, end) when:
        # pos + cn_binsize > start and pos < end.
        left = np.searchsorted(pos, int(start) - cn_binsize, side='right')
        right = np.searchsorted(pos, int(end), side='left')

        if left >= right:
            return np.nan, np.nan

        interval_cn = cn_value[left:right]
        interval_count = count_value[left:right]

        interval_cn = interval_cn[np.isfinite(interval_cn)]
        interval_count = interval_count[np.isfinite(interval_count)]

        median_cn = float(np.median(interval_cn)) if interval_cn.size else np.nan

        mean_count = float(np.mean(interval_count)) if interval_count.size else np.nan

        return median_cn, mean_count

    # ------------------------------------------------------------------
    # Helper functions
    # ------------------------------------------------------------------
    def build_boundary_index(segments, include_segment_index=False):
        index = defaultdict(list)

        for seg_idx, row in segments.iterrows():
            chrom = str(row['Chrom'])

            if include_segment_index:
                index[(chrom, 'L')].append((int(row['Start']), seg_idx))
                index[(chrom, 'R')].append((int(row['End']), seg_idx))
            else:
                index[(chrom, 'L')].append(int(row['Start']))
                index[(chrom, 'R')].append(int(row['End']))

        for key in index:
            index[key] = sorted(index[key])

        return index

    def partner_is_anchored(partner, boundary_index, window, current_seg_idx=None):
        """
        Check whether a PE partner anchors another amplified boundary.
        """
        pchrom, pcoord, pstrand = partner
        pchrom = str(pchrom)
        pcoord = int(pcoord)
        pside = 'R' if pstrand == '+' else 'L'

        candidates = boundary_index.get((pchrom, pside), [])

        for item in candidates:
            if isinstance(item, tuple):
                boundary, seg_idx = item
            else:
                boundary, seg_idx = item, None

            if current_seg_idx is not None and seg_idx == current_seg_idx:
                continue

            if abs(int(boundary) - pcoord) <= window:
                return True

        return False

    def nearest_boundary(chrom, coord, side, window, boundary_index):
        """
        Find the nearest same-direction boundary using a sorted index.
        """
        positions = boundary_index.get((str(chrom), side), [])

        if not positions:
            return None

        coord = int(coord)
        index = np.searchsorted(positions, coord)

        candidates = []

        if index < len(positions):
            candidates.append(int(positions[index]))

        if index > 0:
            candidates.append(int(positions[index - 1]))

        if not candidates:
            return None

        nearest = min(candidates, key=lambda position: abs(position - coord))

        if abs(nearest - coord) > window:
            return None

        return nearest

    # ------------------------------------------------------------------
    # Pass 1: refine boundaries with Clip == 0
    # ------------------------------------------------------------------
    initial_boundary_index = build_boundary_index(out, include_segment_index=True)

    n_refined_zero = 0

    for current_seg_idx, seg in out.iterrows():
        changed = False

        for side, pos_col, clip_col, disc_col, partner_col in [
            ('L', 'Start', 'ClipLeft',  'DiscLeft',  'PartnerLeft'),
            ('R', 'End',   'ClipRight', 'DiscRight', 'PartnerRight'),
        ]:
            clip_count = pd.to_numeric(seg.get(clip_col, np.nan), errors='coerce')

            # Only Clip == 0 is handled in this pass.
            if pd.isnull(clip_count) or float(clip_count) != 0:
                continue

            chrom = str(seg['Chrom'])
            old_coord = int(seg[pos_col])
            candidates = []

            for record in edge_index.get((chrom, side), []):
                (_chrom, coord, _side, support, partner, _edge_id, _end_id) = record

                distance = abs(coord - old_coord)
                if distance > binsize:
                    continue

                if not partner_is_anchored(partner, initial_boundary_index, partner_anchor_window, current_seg_idx):
                    continue

                candidates.append((distance, -support, coord, support, partner))

            if not candidates:
                continue

            _, _, coord, support, partner = min(candidates)

            # Ensure that the adjusted segment remains valid.
            if side == 'L':
                if coord >= int(out.at[current_seg_idx, 'End']):
                    continue
            else:
                if coord <= int(out.at[current_seg_idx, 'Start']):
                    continue

            out.at[current_seg_idx, pos_col] = coord
            out.at[current_seg_idx, disc_col] = support
            out.at[current_seg_idx, partner_col] = partner

            changed = True
            n_refined_zero += 1

        if changed:
            out.at[current_seg_idx, 'SegSource'] = 'PE_boundary'

    # ------------------------------------------------------------------
    # Pass 2: refine positive but weak clip boundaries
    # ------------------------------------------------------------------
    # Use boundaries after the Clip==0 correction as partner anchors.
    weak_boundary_index = build_boundary_index(out, include_segment_index=True)

    n_refined_weak = 0

    for current_seg_idx, seg in out.iterrows():
        left_clip = pd.to_numeric(seg.get('ClipLeft', np.nan), errors='coerce')
        right_clip = pd.to_numeric(seg.get('ClipRight', np.nan), errors='coerce')

        # Both ends must have weak clip evidence.
        if pd.isnull(left_clip) or pd.isnull(right_clip):
            continue

        if float(left_clip) > weak_clip_max or float(right_clip) > weak_clip_max:
            continue

        for side, pos_col, clip_col, disc_col, partner_col in [
            ('L', 'Start', 'ClipLeft',  'DiscLeft',  'PartnerLeft'),
            ('R', 'End',   'ClipRight', 'DiscRight', 'PartnerRight'),
        ]:
            clip_count = pd.to_numeric(out.at[current_seg_idx, clip_col], errors='coerce')

            # Clip==0 was handled above. This pass handles only positive
            # weak clip evidence, normally Clip==1.
            if pd.isnull(clip_count) or float(clip_count) <= 0 or float(clip_count) > weak_clip_max:
                continue

            chrom = str(out.at[current_seg_idx, 'Chrom'])
            old_coord = int(out.at[current_seg_idx, pos_col])
            candidates = []

            for record in edge_index.get((chrom, side), []):
                (_chrom, coord, _side, support, partner, _edge_id, _end_id) = record

                if support < weak_clip_min_edge_support:
                    continue

                distance = abs(coord - old_coord)

                if distance > weak_clip_boundary_window:
                    continue

                # The PE endpoint must extend the segment outward.
                if side == 'L' and coord >= old_coord:
                    continue

                if side == 'R' and coord <= old_coord:
                    continue

                # The other PE endpoint must accurately anchor another
                # amplified segment boundary.
                if not partner_is_anchored(partner, weak_boundary_index, weak_clip_partner_window, current_seg_idx):
                    continue

                if side == 'L':
                    new_start = coord
                    new_end = int(out.at[current_seg_idx, 'End'])
                else:
                    new_start = int(out.at[current_seg_idx, 'Start'])
                    new_end = coord

                new_length = new_end - new_start

                if not (min_segment_size <= new_length <= max_segment_size):
                    continue

                # Check only the newly included interval.
                added_start = min(old_coord, coord)
                added_end = max(old_coord, coord)

                added_cn, _ = interval_stats(chrom, added_start, added_end)

                if not np.isfinite(added_cn) or added_cn < min_cn:
                    continue

                candidates.append((distance, -support, coord, support, partner))

            if not candidates:
                continue

            # Nearest endpoint first; support breaks a distance tie.
            _, _, coord, support, partner = min(candidates)

            out.at[current_seg_idx, pos_col] = coord

            # The old Clip value belongs to the old coordinate.
            out.at[current_seg_idx, clip_col] = 0
            out.at[current_seg_idx, disc_col] = support
            out.at[current_seg_idx, partner_col] = partner
            out.at[current_seg_idx, 'SegSource'] = 'PE_boundary'

            n_refined_weak += 1

    # ------------------------------------------------------------------
    # Build final boundaries after both refinement passes
    # ------------------------------------------------------------------
    final_boundary_index = build_boundary_index(out, include_segment_index=False)

    # ------------------------------------------------------------------
    # Find orphan PE endpoints for missing-segment rescue
    # ------------------------------------------------------------------
    orphan = defaultdict(list)

    for (chrom, coord, side, support, partner, edge_id, end_id) in edge_end:

        # Its own endpoint must be far from every existing boundary.
        if nearest_boundary(chrom, coord, side, orphan_window, final_boundary_index) is not None:
            continue

        pchrom, pcoord, pstrand = partner
        pside = 'R' if pstrand == '+' else 'L'

        # Its partner must anchor an existing amplified boundary.
        if nearest_boundary(pchrom, pcoord, pside, partner_anchor_window, final_boundary_index) is None:
            continue

        orphan[(chrom, side)].append((coord, support, partner, edge_id, end_id))

    def collapse_endpoint(records):
        """
        Collapse nearby PE endpoints and retain the strongest representative.
        """
        records = sorted(records, key=lambda item: item[0])

        collapsed = []
        group = []

        for record in records:
            if group and record[0] - group[-1][0] > edge_merge_window:
                collapsed.append(max(group, key=lambda item: (item[1], -item[0])))
                group = []

            group.append(record)

        if group:
            collapsed.append(max(group, key=lambda item: (item[1], -item[0])))

        return collapsed

    # Existing segment intervals, used to prevent PE_seg overlap.
    existing = defaultdict(list)

    for row in out.itertuples(index=False):
        existing[str(row.Chrom)].append((int(row.Start), int(row.End)))

    def has_overlap(chrom, start, end):
        return any(start < old_end and end > old_start for old_start, old_end in existing.get(str(chrom), []))

    # ------------------------------------------------------------------
    # Rescue missing PE segments
    # ------------------------------------------------------------------
    rescued = []

    orphan_chromosomes = sorted({chrom for chrom, _side in orphan})

    for chrom in orphan_chromosomes:
        left_end = collapse_endpoint(orphan.get((chrom, 'L'), []))
        right_end = collapse_endpoint(orphan.get((chrom, 'R'), []))

        if not left_end or not right_end:
            continue

        right_positions = np.asarray([record[0] for record in right_end], dtype=np.int64)

        used_right = set()

        for left_record in left_end:
            first_right = np.searchsorted(right_positions, left_record[0] + min_segment_size, side='left')

            # Only test the closest few R candidates.
            for right_idx in range(first_right, min(first_right + 5, len(right_end))):

                if right_idx in used_right:
                    continue

                right_record = right_end[right_idx]
                start = int(left_record[0])
                end = int(right_record[0])
                length = end - start

                if length > max_segment_size:
                    break

                if has_overlap(chrom, start, end):
                    continue

                median_cn, mean_count = interval_stats(chrom, start, end)

                if not np.isfinite(median_cn) or median_cn < min_cn:
                    continue

                rescued.append((chrom, left_record, right_record, median_cn, mean_count))

                existing[chrom].append((start, end))
                used_right.add(right_idx)
                break

    # ------------------------------------------------------------------
    # Add rescued PE segments
    # ------------------------------------------------------------------
    if rescued:
        rescued_rows = []

        for (chrom, left_record, right_record, median_cn, mean_count) in rescued:

            row = {column: np.nan for column in out.columns}

            row.update({
                'Chrom': chrom,
                'Start': int(left_record[0]),
                'End': int(right_record[0]),
                'CN': median_cn,
                'ClipLeft': 0,
                'ClipRight': 0,
                'DiscLeft': int(left_record[1]),
                'DiscRight': int(right_record[1]),
                'PartnerLeft': left_record[2],
                'PartnerRight': right_record[2],
                'SegSource': 'PE_seg',
            })

            if 'Count' in row:
                row['Count'] = mean_count

            if 'Gene' in row:
                row['Gene'] = ''

            if 'CancerGene' in row:
                row['CancerGene'] = ''

            rescued_rows.append(row)

        out = pd.concat([out, pd.DataFrame(rescued_rows)], ignore_index=True)

    out = out.sort_values(['Chrom', 'Start', 'End']).reset_index(drop=True)

    # ------------------------------------------------------------------
    # Snap PE endpoints to final segment boundaries
    # ------------------------------------------------------------------
    snap_index = build_boundary_index(out, include_segment_index=False)

    edges_new = []
    n_snapped = 0

    for edge in edges:
        if len(edge) < 3:
            edges_new.append(edge)
            continue

        new_endpoints = []

        for chrom, coord, strand in (edge[0], edge[1]):
            side = 'R' if strand == '+' else 'L'

            target = nearest_boundary(chrom, coord, side, edge_merge_window, snap_index)

            if target is not None and target != int(coord):
                coord = target
                n_snapped += 1

            new_endpoints.append((chrom, int(coord), strand))

        # READID: boundary snapping changes coordinates only; keep support IDs/barcodes.
        edges_new.append((new_endpoints[0], new_endpoints[1], *edge[2:]))

    print(
        f'refine zero-clip boundary: {n_refined_zero}; '
        f'refine weak-clip boundary: {n_refined_weak}; '
        f'rescue missing PE segment: {len(rescued)}; '
        f'snap PE endpoints: {n_snapped};'
    )

    return out, edges_new


#------------------------------------------------------------------------------#
# refine amplified segments with discordant PE endpoints
#------------------------------------------------------------------------------#
def refine_segment_with_discordant(
        cn_amp, bp_all, edges, bin_count,
        bp_min_clip = 2, pre_window = 150, boundary_window = 10000,
        min_edge_support = 2, min_cn = 4.5,
        min_segment_size = 100, max_segment_size = 100000000,
        binsize = 10000):
    """Use discordant PE endpoints to refine amplified segments.

    Existing segment:
        always keep the original clip-refined segment. If clip support on a
        boundary is below bp_min_clip, add one PE-refined alternative segment
        instead of overwriting the original segment.
    Missing segment:
        pair orphan L/R PE endpoints on the same chromosome, require both PE
        partners to anchor existing amplified boundaries, and require median
        unsmoothed bin CN >= min_cn before adding the interval.

    The returned edges are snapped to repaired segment boundaries within
    pre_window, so downstream (Chrom, Coord, Clip) exact graph matching uses
    one canonical coordinate.

    Returns
    -------
    cn_amp_new, edges_new
    """
    out = cn_amp.copy().reset_index(drop = True)
    edges = list(edges) if edges is not None else []

    for col in ['DiscLeft', 'DiscRight', 'PartnerLeft', 'PartnerRight']:
        if col not in out.columns:
            out[col] = 0 if 'Disc' in col else None
    if 'SegSource' not in out.columns:
        out['SegSource'] = 'clip'
    else:
        out['SegSource'] = out['SegSource'].fillna('clip')

    # Build a sorted endpoint index once: O(E log E). '+' is R and '-' is L.
    edge_index = defaultdict(list)
    edge_end = []
    for edge_id, edge in enumerate(edges):
        if len(edge) < 3:
            continue
        A, B, n = edge[0], edge[1], int(edge[2])
        if n < min_edge_support:
            continue
        for end_id, (this_end, other_end) in enumerate(((A, B), (B, A))):
            chrom, coord, strand = this_end
            clip = 'R' if strand == '+' else 'L'
            rec = (int(coord), n, other_end, edge_id, end_id)
            edge_index[(str(chrom), clip)].append(rec)
            edge_end.append((str(chrom), int(coord), clip, n,
                             other_end, edge_id, end_id))

    for key in edge_index:
        edge_index[key].sort(key = lambda x: x[0])

    edge_coord = {
        key: np.asarray([x[0] for x in value], dtype = np.int64)
        for key, value in edge_index.items()
    }

    def nearest_edge(chrom, coord, clip, window):
        """Strongest nearby edge; distance is the tie breaker."""
        key = (str(chrom), clip)
        pos = edge_coord.get(key)
        if pos is None or pos.size == 0:
            return None
        left = np.searchsorted(pos, int(coord) - window, side = 'left')
        right = np.searchsorted(pos, int(coord) + window, side = 'right')
        if left == right:
            return None

        cand = edge_index[key][left:right]
        return max(cand, key = lambda x: (x[1], -abs(x[0] - int(coord))))

    # Build the clip index before refining existing segments because the
    # refinement loop uses nearby_clip_count() to annotate PE alternatives.
    clip_index = defaultdict(list)
    if bp_all is not None:
        for side in ('left', 'right'):
            df = bp_all.get(side)
            if df is None or df.empty:
                continue
            for row in df.itertuples(index = False):
                clip_index[(str(row.Chrom), str(row.Clip))].append(
                    (int(row.Coord), int(row.Count)))
    for key in clip_index:
        clip_index[key].sort(key = lambda x: x[0])

    clip_coord = {
        key: np.asarray([x[0] for x in value], dtype = np.int64)
        for key, value in clip_index.items()
    }

    def nearby_clip_count(chrom, coord, clip):
        key = (str(chrom), clip)
        pos = clip_coord.get(key)
        if pos is None or pos.size == 0:
            return 0
        left = np.searchsorted(pos, int(coord) - pre_window, side = 'left')
        right = np.searchsorted(pos, int(coord) + pre_window, side = 'right')
        if left == right:
            return 0
        return max(x[1] for x in clip_index[key][left:right])

    # Keep every clip-refined segment and add at most one PE alternative per row.
    n_refined = 0
    pe_segment = []
    for _, seg in out.iterrows():
        pe_seg = seg.copy()
        changed = False
        for side, pos_col, clip_col, disc_col, partner_col in [
            ('L', 'Start', 'ClipLeft', 'DiscLeft', 'PartnerLeft'),
            ('R', 'End',   'ClipRight', 'DiscRight', 'PartnerRight'),
        ]:
            clip_count = pd.to_numeric(seg.get(clip_col, 0), errors = 'coerce')
            clip_count = 0 if pd.isnull(clip_count) else float(clip_count)
            if clip_count >= bp_min_clip:
                continue

            best = nearest_edge(seg['Chrom'], seg[pos_col], side,
                                boundary_window)
            if best is None:
                continue
            coord, support, partner, _, _ = best
            if side == 'L' and coord >= int(seg['End']):
                continue
            if side == 'R' and coord <= int(seg['Start']):
                continue
            pe_seg[pos_col] = coord
            pe_seg[clip_col] = nearby_clip_count(seg['Chrom'], coord, side )
            pe_seg[disc_col] = support
            pe_seg[partner_col] = partner
            changed = True
            n_refined += 1

        if changed and (int(pe_seg['Start']) != int(seg['Start']) or
                        int(pe_seg['End']) != int(seg['End'])):
            pe_seg['SegSource'] = 'PE_boundary'
            pe_segment.append(pe_seg)

    if pe_segment:
        out = pd.concat([out, pd.DataFrame(pe_segment)], ignore_index = True)

    def build_boundary_index(df):
        index = defaultdict(list)
        for row in df.itertuples(index = False):
            index[(str(row.Chrom), 'L')].append(int(row.Start))
            index[(str(row.Chrom), 'R')].append(int(row.End))
        return {key: np.sort(np.asarray(value, dtype = np.int64))
                for key, value in index.items()}

    boundary_index = build_boundary_index(out)

    def nearest_boundary(chrom, coord, clip, window, index = None):
        index = boundary_index if index is None else index
        pos = index.get((str(chrom), clip))
        if pos is None or pos.size == 0:
            return None
        i = np.searchsorted(pos, int(coord))
        cand = []
        if i < pos.size:
            cand.append(int(pos[i]))
        if i > 0:
            cand.append(int(pos[i - 1]))
        if not cand:
            return None
        best = min(cand, key = lambda x: abs(x - int(coord)))
        return best if abs(best - int(coord)) <= window else None

    # Use unsmoothed bin CN for missing-segment validation.
    cn_index = {}
    if bin_count is not None and not bin_count.empty:
        for chrom, df in bin_count.groupby('Chrom', sort = False):
            df = df.sort_values('Coord')
            cn_index[str(chrom)] = (
                df['Coord'].to_numpy(dtype = np.int64),
                pd.to_numeric(df['CN'], errors = 'coerce').to_numpy(dtype = float),
                pd.to_numeric(df['Count'], errors = 'coerce').to_numpy(dtype = float))

    def interval_stats(chrom, start, end):
        data = cn_index.get(str(chrom))
        if data is None:
            return np.nan, np.nan
        pos, cn_value, count_value = data
        left = np.searchsorted(pos, int(start) - binsize, side = 'left')
        right = np.searchsorted(pos, int(end), side = 'right')
        cn_value = cn_value[left:right]
        count_value = count_value[left:right]
        cn_value = cn_value[np.isfinite(cn_value)]
        count_value = count_value[np.isfinite(count_value)]
        cn = float(np.median(cn_value)) if cn_value.size > 0 else np.nan
        count = float(np.mean(count_value)) if count_value.size > 0 else np.nan
        return cn, count

    # Orphan endpoint = its own side is not represented by cn_amp, while the
    # other end of the same PE edge is anchored to an existing segment boundary.
    orphan = defaultdict(list)
    for chrom, coord, clip, support, partner, edge_id, end_id in edge_end:
        if nearest_boundary(chrom, coord, clip, boundary_window) is not None:
            continue
        pchrom, pcoord, pstrand = partner
        pclip = 'R' if pstrand == '+' else 'L'
        if nearest_boundary(pchrom, pcoord, pclip, boundary_window) is None:
            continue
        orphan[(chrom, clip)].append(
            (coord, support, partner, edge_id, end_id))

    # Collapse endpoints from nearby edge clusters; retain the strongest one.
    def collapse_endpoint(value):
        value = sorted(value, key = lambda x: x[0])
        op, group = [], []
        for rec in value:
            if group and rec[0] - group[-1][0] > pre_window:
                op.append(max(group, key = lambda x: (x[1], -x[0])))
                group = []
            group.append(rec)
        if group:
            op.append(max(group, key = lambda x: (x[1], -x[0])))
        return op

    existing = defaultdict(list)
    for row in out.itertuples(index = False):
        existing[str(row.Chrom)].append((int(row.Start), int(row.End)))

    def has_overlap(chrom, start, end):
        return any(start < e and end > s for s, e in existing.get(str(chrom), []))

    rescued = []
    for chrom in sorted(set(k[0] for k in orphan)):
        left_end = collapse_endpoint(orphan.get((chrom, 'L'), []))
        right_end = collapse_endpoint(orphan.get((chrom, 'R'), []))
        if not left_end or not right_end:
            continue
        right_pos = np.asarray([x[0] for x in right_end], dtype = np.int64)
        used_right = set()

        for L in left_end:
            i = np.searchsorted(right_pos, L[0] + min_segment_size,
                                side = 'left')
            # Check only the nearest few right endpoints to prevent L x R blow-up.
            for j in range(i, min(i + 5, len(right_end))):
                R = right_end[j]
                if j in used_right or R[0] - L[0] > max_segment_size:
                    continue
                if has_overlap(chrom, L[0], R[0]):
                    continue
                cn, count = interval_stats(chrom, L[0], R[0])
                if not np.isfinite(cn) or cn < min_cn:
                    continue
                rescued.append((chrom, L, R, cn, count))
                existing[chrom].append((L[0], R[0]))
                used_right.add(j)
                break

    if rescued:
        rows = []
        for chrom, L, R, cn, count in rescued:
            row = {col: np.nan for col in out.columns}
            row.update({
                'Chrom': chrom, 'Start': L[0], 'End': R[0], 'CN': cn,
                'ClipLeft': nearby_clip_count(chrom, L[0], 'L'),
                'ClipRight': nearby_clip_count(chrom, R[0], 'R'),
                'DiscLeft': L[1], 'DiscRight': R[1],
                'PartnerLeft': L[2], 'PartnerRight': R[2],
                'SegSource': 'PE_seg',
            })
            if 'Count' in row:
                row['Count'] = count
            if 'Gene' in row:
                row['Gene'] = ''
            if 'CancerGene' in row:
                row['CancerGene'] = ''
            rows.append(row)
        out = pd.concat([out, pd.DataFrame(rows)], ignore_index = True)

    out = out.sort_values(['Chrom', 'Start', 'End']).reset_index(drop = True)

    # Canonicalize PE coordinates to the final segment boundary within 150 bp.
    final_boundary_index = build_boundary_index(out)
    edges_new, n_snapped = [], 0
    for edge in edges:
        if len(edge) < 3:
            edges_new.append(edge)
            continue
        A, B = edge[0], edge[1]
        new_end = []
        for chrom, coord, strand in (A, B):
            clip = 'R' if strand == '+' else 'L'
            target = nearest_boundary(chrom, coord, clip, pre_window,
                                      index = final_boundary_index)
            if target is not None and target != int(coord):
                coord = target
                n_snapped += 1
            new_end.append((chrom, int(coord), strand))
        # READID: segment refinement must not drop compact support evidence.
        edges_new.append((new_end[0], new_end[1], *edge[2:]))

    print(f'keep clip segment; add PE boundary segment: {len(pe_segment)} '
          f'({n_refined} boundaries); '
          f'rescue missing segment: {len(rescued)}; '
          f'snap PE endpoint: {n_snapped}')
    return out, edges_new

#------------------------------------------------------------------------------#
# refine discordant PE endpoints with SA reads
#------------------------------------------------------------------------------#

import numpy as np
import pandas as pd
from collections import defaultdict


def _clip_of(strand):
    if strand == '+':
        return 'R'
    if strand == '-':
        return 'L'
    raise ValueError(f"invalid strand value: {strand!r}, expected '+' or '-'")


def _build_pair_index(duo):
    """Build key -> (coordsA sorted, coordsB, support) for fast lookup.

    key = (ChromA, ClipA, ChromB, ClipB). Each row is stored twice
    (forward and reversed) so edge endpoint order does not matter.
    Arrays are sorted by CoordA to allow searchsorted windowing.
    """
    raw = defaultdict(lambda: ([], [], []))

    for c1, p1, x1, c2, p2, x2, s in zip(
            duo['Chrom1'], duo['Clip1'], duo['Coord1'],
            duo['Chrom2'], duo['Clip2'], duo['Coord2'],
            duo['Count']):
        a, b, sup = raw[(c1, p1, c2, p2)]
        a.append(x1); b.append(x2); sup.append(s)
        a, b, sup = raw[(c2, p2, c1, p1)]
        a.append(x2); b.append(x1); sup.append(s)

    index = {}
    for key, (a, b, sup) in raw.items():
        a = np.asarray(a, dtype=np.int64)
        order = np.argsort(a, kind='stable')
        index[key] = (
            a[order],
            np.asarray(b, dtype=np.int64)[order],
            np.asarray(sup, dtype=float)[order],
        )
    return index


def _best_match(bucket, coord_a, coord_b, max_distance):
    """Find the highest-support candidate within max_distance of both ends.

    Returns (new_coord_a, new_coord_b), 'ambiguous' (tied top support),
    or None (no candidate in range).
    """
    coords_a, coords_b, supports = bucket

    lo = np.searchsorted(coords_a, coord_a - max_distance, side='left')
    hi = np.searchsorted(coords_a, coord_a + max_distance, side='right')
    if lo == hi:
        return None

    mask = np.abs(coords_b[lo:hi] - coord_b) <= max_distance
    if not mask.any():
        return None

    cand_a = coords_a[lo:hi][mask]
    cand_b = coords_b[lo:hi][mask]
    cand_s = supports[lo:hi][mask]

    # Same exact coordinate pair may appear more than once; keep max support.
    best_by_coord = {}
    for a, b, s in zip(cand_a, cand_b, cand_s):
        k = (int(a), int(b))
        if k not in best_by_coord or s > best_by_coord[k]:
            best_by_coord[k] = s

    top = max(best_by_coord.values())
    winners = [k for k, s in best_by_coord.items() if s == top]

    return winners[0] if len(winners) == 1 else 'ambiguous'


def _edge_read_id_array(value):
    """Normalize compact ReadID arrays and upgrade legacy QNAME edge payloads."""
    if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.integer):
        return value.astype(np.uint64, copy=False)

    values = list(value) if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)) else []
    if not values:
        return np.empty(0, dtype=np.uint64)
    if all(isinstance(item, (int, np.integer)) for item in values):
        return np.asarray(values, dtype=np.uint64)

    # READID: compatibility for old edges that still carry QNAME strings.
    read_ids, _ = covbam.encode_read_evidence(values)
    return read_ids


def _collapse_by_endpoints(refined_edges):
    collapsed = {}
    has_read_ids = any(len(edge) > 3 for edge in refined_edges)
    has_readbarcodes = any(len(edge) > 4 for edge in refined_edges)

    for edge in refined_edges:
        if len(edge) < 3:
            continue

        end_a, end_b, support = edge[:3]
        read_ids = _edge_read_id_array(edge[3] if len(edge) > 3 else [])
        readbarcodes = edge[4] if len(edge) > 4 else []
        readbarcodes = list(readbarcodes) if isinstance(readbarcodes, (list, tuple, np.ndarray, pd.Series)) else []

        key = (end_a, end_b) if end_a <= end_b else (end_b, end_a)
        if key not in collapsed:
            collapsed[key] = {
                'end_a': key[0],
                'end_b': key[1],
                'read_ids': [],
                'readbarcodes': [],
                'fallback_support': 0,
            }

        for index, read_id in enumerate(read_ids):
            collapsed[key]['read_ids'].append(read_id)
            barcode = readbarcodes[index] if index < len(readbarcodes) else None
            barcode = str(barcode) if barcode is not None and str(barcode) not in {'', 'None', 'nan'} else None
            collapsed[key]['readbarcodes'].append(barcode)

        collapsed[key]['fallback_support'] += int(support)

    edges_new = []
    for item in collapsed.values():
        read_ids = np.asarray(item['read_ids'], dtype=np.uint64)
        readbarcodes = item['readbarcodes']
        if not any(barcode is not None for barcode in readbarcodes):
            readbarcodes = []
        support = item['fallback_support']

        # BARCODE: collapsed WGS edges retain no per-read None placeholders.
        if has_readbarcodes:
            edges_new.append((item['end_a'], item['end_b'], support, read_ids, readbarcodes))
        elif has_read_ids:
            edges_new.append((item['end_a'], item['end_b'], support, read_ids))
        else:
            edges_new.append((item['end_a'], item['end_b'], support))

    edges_new.sort(key=lambda edge: -edge[2])
    return edges_new, len(refined_edges) - len(edges_new)


def snap_discordant_edges_to_bp_duo(
        edges,
        bp_duo_bam,
        max_distance=150,
        min_support=1):
    """Refine discordant PE edges using split-read breakpoint pairs.

    Rules
    -----
    1. strand '+' -> Clip 'R', strand '-' -> Clip 'L' (anything else raises).
    2. Both endpoints of a PE edge must match the same breakpoint pair,
       each within max_distance.
    3. Pick the candidate with the highest support; ties keep the
       original edge unchanged.
    4. Edges colliding onto the same endpoint pair after refinement merge
       their ReadID occurrences and aligned cell barcodes; support counts are
       added without read-level deduplication.
    """

    edges = list(edges) if edges is not None else []

    empty_stats = {
        'snapped_endpoints': 0,
        'snapped_edges': 0,
        'merged_edges': 0,
        'ambiguous_edges': 0,
        'no_match_edges': len(edges),
    }

    if not edges or bp_duo_bam is None or bp_duo_bam.empty:
        return edges, empty_stats

    required = {'Chrom1', 'Coord1', 'Clip1', 'Chrom2', 'Coord2', 'Clip2', 'Count'}
    missing = required.difference(bp_duo_bam.columns)
    if missing:
        raise ValueError(f"bp_duo_bam missing required columns: {sorted(missing)}")

    duo = bp_duo_bam.loc[:, list(required)].copy()

    for col in ('Chrom1', 'Chrom2', 'Clip1', 'Clip2'):
        duo[col] = duo[col].astype(str)
    duo['Coord1'] = pd.to_numeric(duo['Coord1'], errors='coerce')
    duo['Coord2'] = pd.to_numeric(duo['Coord2'], errors='coerce')
    duo['Count'] = pd.to_numeric(duo['Count'], errors='coerce').fillna(0)

    duo = duo[
        duo['Coord1'].notna() & duo['Coord2'].notna()
        & duo['Clip1'].isin(['L', 'R']) & duo['Clip2'].isin(['L', 'R'])
    ].reset_index(drop=True)

    if duo.empty:
        return edges, empty_stats

    duo['Coord1'] = duo['Coord1'].astype(np.int64)
    duo['Coord2'] = duo['Coord2'].astype(np.int64)

    duo = duo[duo['Count'] >= min_support]
    if duo.empty:
        return edges, empty_stats

    pair_index = _build_pair_index(duo)

    refined = []
    n_snapped_endpoints = n_snapped_edges = 0
    n_ambiguous = n_no_match = 0

    for edge in edges:
        if len(edge) < 3:
            refined.append(edge)
            continue

        end_a, end_b, pe_support = edge[:3]
        chrom_a, coord_a, strand_a = end_a
        chrom_b, coord_b, strand_b = end_b
        coord_a, coord_b = int(coord_a), int(coord_b)

        key = (str(chrom_a), _clip_of(strand_a), str(chrom_b), _clip_of(strand_b))
        bucket = pair_index.get(key)
        result = _best_match(bucket, coord_a, coord_b, max_distance) if bucket else None

        if result is None:
            refined.append(edge)
            n_no_match += 1
            continue

        if result == 'ambiguous':
            refined.append(edge)
            n_ambiguous += 1
            continue

        new_coord_a, new_coord_b = result
        changed_a = new_coord_a != coord_a
        changed_b = new_coord_b != coord_b
        n_snapped_endpoints += int(changed_a) + int(changed_b)
        n_snapped_edges += int(changed_a or changed_b)

        # READID: SA snapping preserves compact IDs for later union-by-junction.
        refined.append(((chrom_a, new_coord_a, strand_a), (chrom_b, new_coord_b, strand_b), *edge[2:]))

    edges_new, n_merged = _collapse_by_endpoints(refined)

    stats = {
        'snapped_endpoints': n_snapped_endpoints,
        'snapped_edges': n_snapped_edges,
        'merged_edges': n_merged,
        'ambiguous_edges': n_ambiguous,
        'no_match_edges': n_no_match,
    }
    return edges_new, stats
