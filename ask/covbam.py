################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd
import math
import logging
import os
import time
import random
from collections import defaultdict

#------------------------------------------------------------------------------#
from grange import GRange
import misc

#------------------------------------------------------------------------------#
def _alignment_filename(alignment_file):
    """Return a pysam AlignmentFile filename as text."""
    filename = alignment_file.filename
    if isinstance(filename, bytes):
        filename = filename.decode()
    return filename

#------------------------------------------------------------------------------#
def downsampled_bam_path(bamfile, downsample_ratio):
    """Build the output path for a downsampled BAM."""
    bam_pathname = _alignment_filename(bamfile.bamfile)
    root, ext = os.path.splitext(bam_pathname)
    ratio_tag = f'{downsample_ratio:.2f}'.rstrip('0').rstrip('.')

    if ext.lower() == '.bam':
        return f'{root}_{ratio_tag}_downsample.bam'
    return f'{bam_pathname}_{ratio_tag}_downsample.bam'

#------------------------------------------------------------------------------#
# function to check read
#------------------------------------------------------------------------------#
def check_read(read, mapq = 20, nmmax = 1):
    """
    true if the read meets certian conditions
    """
    return not read.is_unmapped \
        and read.mapping_quality >= mapq \
        and read.get_tag('NM') <= nmmax \
        and not read.is_duplicate

#------------------------------------------------------------------------------#
# get reads  with SA from bamfile
#------------------------------------------------------------------------------#
def get_reads_with_SA(read, reads_with_SA, nmmax=1,  SA_with_nm = True):  
    '''
    get bp pair reads from bamfile
    '''
    tf = True
    if SA_with_nm: 
        tf = check_read(read, mapq = 0, nmmax = nmmax)   
    if read.has_tag('SA') and tf:
        # read name
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
            if 'single' not in reads_with_SA[read_name]:
                reads_with_SA[read_name]['single'] = []
            reads_with_SA[read_name]['single'].append(read)
    return reads_with_SA

#------------------------------------------------------------------------------#
# get all clip depth
#------------------------------------------------------------------------------#
def clip_from_bam(bamfile, mapq = 20, nmmax = 1, n = 0,  SA_with_nm = True):
    """
    get all clip positions from bam files
    -- count softclip and hardclip reads per chromosome coordinates
    for left and right clip, respectively

    usage:
        clip_from_bam(bamfile)

    input: bamfile
    output: pandas dataframe of split read counts in a dict
    nmmax: max number of mismatches
    mapq: minimal mapping quality
    readq: minimal read mean quality (readq = 20, not necessary)
    """

    # with pysam.AlignmentFile(bamfile, "rb") as bamf:

    left = list() # init left clip
    right = list() # init right clip
    reads_with_SA = {}  # init reads with SA

    # save downsampled reads for downstream breakpoint-pair analysis
    outfile = None
    if bamfile.downsample_ratio != 1:
        ds_bamfile_pathname = downsampled_bam_path(bamfile, bamfile.downsample_ratio)
        print(f'writing downsampled BAM: {ds_bamfile_pathname}')
        outfile = pysam.AlignmentFile(ds_bamfile_pathname, 'wb', template=bamfile.bamfile)

    try:
        for read in bamfile.fetch():
            if outfile is not None:
                outfile.write(read)
            reads_with_SA = get_reads_with_SA(read, reads_with_SA, nmmax, SA_with_nm)
            if check_read(read, mapq = mapq, nmmax = nmmax):

                # left clip read
                if (read.cigartuples[0][0] == 4): # soft clip reads
                    left.append((read.reference_name, \
                                    read.reference_start, "L"))
                elif (read.cigartuples[0][0] == 5): # hard clip reads
                    left.append((read.reference_name, \
                                    read.reference_start, "L"))

                # right clip reads
                if (read.cigartuples[-1][0] == 4): # soft clip reads
                    right.append((read.reference_name, \
                                    read.reference_end-1, "R"))
                elif (read.cigartuples[-1][0] == 5): # hard clip reads
                    right.append((read.reference_name, \
                                    read.reference_end-1, "R"))
    finally:
        if outfile is not None:
            outfile.close()

    # convert all clip read counts to df
    left = list2df_count(left, n = n, f_sort = False)
    right = list2df_count(right, n = n, f_sort = False)

    # return left and right clip points by chrom.position
    return dict(left=left, right=right), reads_with_SA

#------------------------------------------------------------------------------#
def list2df_count(lst, colnames=["Chrom", "Coord", "Clip"], n = 10, f_sort=True):
    """
    convert list to dataframe, count and sort

    usage:
        list2df_count(lst)

    input: list of softclip reads with "Chrom", "Coord" and "Clip"
    output: pandas dataframe of split read counts

    """
    df = pd.DataFrame(lst, columns=colnames)
    df = df.groupby(colnames).size().reset_index(name='Count')
    if (n > 0):
        df = df[df['Count']>=n]
    if (f_sort):
        df = df.sort_values('Count', ascending=False)
    return df

#------------------------------------------------------------------------------#
def clip2bedgraph(clip, n = 2):
    """
    clip to bedgraph for visualization
    """
    # init
    left = clip['left'].copy()
    right = clip['right'].copy()

    # filter
    left = left[left['Count'] >= n]
    right = right[right['Count'] >= n]

    # output bedgraph format to a pandas dataframe
    left['Start'] = left['Coord']
    left['End'] = left['Start'] + 1
    right['Start'] = right['Coord']
    right['End'] = right['Start'] + 1
    left['Count'] = -left['Count']
    df = pd.concat([left[['Chrom', 'Start', 'End', 'Count']],
                    right[['Chrom', 'Start', 'End', 'Count']]])
    return df.sort_values(['Chrom', 'Start', 'End'])

# DFS  ################# ---------------------------------------

#------------------------------------------------------------------------------#
# apply heavy smoothing on sub bin counts to remove peaks in ATAC-seq
#------------------------------------------------------------------------------#
def smooth_count(bin_count, k = 20, outlier_sd = 2, npass = 10):
    """smooth bin counts
    """
    df = bin_count.copy()

    df_smooth = []
    for chr in misc.unique(df['Chrom']):
        dfsub = df[df['Chrom'] == chr]
        x_smooth = smooth_count_1d(dfsub['Count'], k, outlier_sd, npass)
        dfsub = dfsub.assign(Count = x_smooth)
        df_smooth.append(dfsub)
    if df_smooth:
        return pd.concat(df_smooth)
    else:
        return pd.DataFrame()


#------------------------------------------------------------------------------#
def smooth_count_1d(x, k = 20, outlier_sd = 2, npass = 10, edge_repeat = False):
    """smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.array(x, dtype = 'float64')
    n = 0

    while n < npass:
        # contruct working matrix
        y = np.zeros((len(x), 2*k + 1), dtype = x.dtype)
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

        # identify outlier
        ym = y[:,k]
        ys = y[:,[i for i in range(2*k+1) if i != k]]
        ys_mean = np.nanmean(ys, axis = 1)
        ys_sd = np.nanstd(ys, ddof=1, axis = 1)
        ys_upper = ys_mean + outlier_sd * ys_sd
        outlier = (ym > ys_upper)

        # assign outlier to na
        x[outlier] = np.nan

        # iteration counter
        n = n + 1

    return x

#------------------------------------------------------------------------------#
# count reads in genomic bins
#------------------------------------------------------------------------------#
def region_count(bamfile, genomesizefile = None,
                 binsize = 10000, mapq = 20,
                 nmmax = 1, sort_ = True):
    """
    count reads in genome bins
    """
    cnt = list()

    # with pysam.AlignmentFile(bamfile, "rb") as bamf:
    for read in bamfile.fetch():
        if check_read(read, mapq = mapq, nmmax = nmmax):
            pos = math.floor(read.reference_start/binsize)*binsize
            cnt.append([read.reference_name, pos])

    # convert all read counts to df
    colnames = ['Chrom', 'Coord']
    df = pd.DataFrame(cnt, columns=colnames)
    df = df.groupby(colnames).size().reset_index(name='Count')
    df['CN'] = 2*df['Count']/np.mean(df['Count'])

    if (sort_): 
        df = df.sort_values(colnames)
    return df


#------------------------------------------------------------------------------#
def region_count_2pass(bamfile, genomesizefile = None,
                   binsize = 10000, mapq = 20,
                   nmmax = 1, sort_ = True,
                   sub_binsize = 100, q = 0.5,
                   colnames = ['Chrom', 'Coord']):
    """wrapper of calculating robust read counts
    in target genomic bins

    Parameters
    ----------
    binsize : target genomic bin size to calculate read counts
    sub_binsize : smaller bin size within target bin size to
        calculate robust score
    q : quantile as robust score, if None auto determine
    """
    # calculate read counts in smaller bin size
    df_sub = region_count(bamfile, genomesizefile, sub_binsize,
                          mapq, nmmax, False)

    # smooth the sub bin counts to remove spikes
    df_sub_smooth = smooth_count(df_sub, k = 2*int(binsize/sub_binsize))

    # calculate read counts in target bin size
    df = region_count_bias(df_sub_smooth, binsize, sub_binsize, q, colnames)

    # fill NA with linear average of nearby counts
    df = df.assign(
        Count = df.groupby(['Chrom'])['Count'].transform(
            lambda group: group.interpolate()
        )
    )

    # calculate the absolute copy number in each bin
    df['CN'] = 2*df['Count']/np.mean(df['Count'])

    if (sort_): 
        df = df.sort_values(colnames)
    return df

#------------------------------------------------------------------------------#
def region_count_bias(df_sub, binsize = 10000, sub_binsize = 100,
                  q = 0.5, colnames = ['Chrom', 'Coord']):
    """calculate robust read counts in target genomic bins
    by using quantiles of read counts in smaller bins
    """
    df = df_sub.copy()
    df = df.assign(Coord = (
        np.floor(np.array(df['Coord'])/binsize)*binsize).astype(int))
    df = df.groupby(colnames)['Count'].quantile(q = q).reset_index(name='Count')
    return df



start = time.time()

# length of chr
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))
DATA_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))

_chrlen_candidates = [
    os.path.join(DATA_REPO, 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai'),
    os.path.join(DATA_REPO, 'hg38.genome'),
    os.path.join(DATA_DIR, 'hg38.genome'),
]
_centromere_candidates = [
    os.path.join(DATA_REPO, 'GRCh38_centromere.bed'),
    os.path.join(DATA_DIR, 'GRCh38_centromere.bed'),
]
_conserved_region_candidates = [
    os.path.join(DATA_REPO, 'conserved_gain5_hg38.bed'),
    os.path.join(DATA_DIR, 'conserved_gain5_hg38.bed'),
]

def _first_existing_file(candidates):
    for filename in candidates:
        if filename and os.path.exists(filename):
            return filename
    return None

def _load_optional_bed(filename, label):
    if filename is None:
        logging.warning("#TIME " + '%.3f\t'%(time.time() - start) + f"No {label} BED file found; using an empty region set.")
        return GRange([], 'list')
    try:
        return GRange(filename, 'bedfile')
    except OSError:
        logging.warning("#TIME " + '%.3f\t'%(time.time() - start) + f"Unable to open {label} BED file: \"{filename}\"")
        return GRange([], 'list')

chrlen_file = _first_existing_file(_chrlen_candidates)
centromere_filename = _first_existing_file(_centromere_candidates)
conserved_regions_filename = _first_existing_file(_conserved_region_candidates)

conserved_regions = _load_optional_bed(conserved_regions_filename, 'conserved regions')
centromere_list = _load_optional_bed(centromere_filename, 'centromere')

# Handling chromosome names, lengths, sorting, positions and addition of new chromosomes
chr_id = {}
chrName = {}
chromList = [str(x) for x in range(1, 23)] + ['X', 'Y']  # must be updated if including an organism with more chroms.

# get chr id
def chrNum(chrname, mode='append'):  # use
    if chrname in chr_id:
        return chr_id[chrname] 
    else:
        if mode == 'init':
            cnum = len(chr_id)
        else:
            cnum = 1000000 + len(chr_id)
        chr_id[chrname] = cnum
        chrName[cnum] = chrname
        return chr_id[chrname]

chrLen = defaultdict(lambda: 0, {})

if chrlen_file is None:
    logging.warning("#TIME " + '%.3f\t'%(time.time() - start) + "No chromosome lengths file found.")
else:
    try:
        with open(chrlen_file) as infile:
            for line in infile:
                ll = line.strip().split()
                if len(ll) >= 2:
                    chrLen[chrNum(ll[0], mode='init')] = int(ll[1])
    except OSError:
        logging.warning("#TIME " + '%.3f\t'%(time.time() - start) + " Unable to open chromosome lengths file: \"" + chrlen_file + "\"")

# get the absolute position given a position on a chr
chrOffset = {} 
def absPos(chrname, pos=0): # use 
    cnum = chrNum(chrname)  # chr7: 6
    if chrNum(chrname) not in chrOffset:  # chrOffset: absolute poisition of start
        chrkeys = sorted(chrName.keys())
        sumlen = sum([chrLen[c] for c in chrLen if c in chrOffset])

        for i in range(len(chrkeys)):
            if chrkeys[i] not in chrOffset:
                chrOffset[chrkeys[i]] = sumlen
                sumlen += chrLen[chrkeys[i]]
            if cnum < chrkeys[i]:
                break
    return chrOffset[chrNum(chrname)] + pos 

def chrPos(abspos): 
    for c in chrOffset: 
        if chrOffset[c] < abspos and chrOffset[c] + chrLen[c] >= abspos:
            return (chrName[c], abspos - chrOffset[c])  
    return None

# def update_chrLen(len_list):
#     for l in len_list:
#         chrLen[chrNum(l[0])] = int(l[1])

# for c in chrLen:
#     ap = absPos(chrName[c]) 

## reference from AA
class prebamfile():
    def __init__(self, bamfile, window_size=10000, read_length=100, max_insert=400, 
                 insert_size=300, min_coverage=30, downsample=10, num_sdevs=3, seed=1):
        self.bamfile = bamfile
        self.window_size = window_size
        self.read_length = read_length
        self.max_insert = max_insert
        self.insert_size = insert_size
        self.min_coverage = min_coverage
        self.downsample = downsample
        self.basic_stats_set = None
        self.seed = seed
        self.interval_coverage_calls = {}
        self.mapping_quality_cutoff = 5
        self.downsample_ratio = 1
        self.num_sdevs = num_sdevs 
        bamfile_pathname = str(self.bamfile.filename.decode()) 
        self.bamfile_filesize = os.path.getsize(bamfile_pathname) /1024**3 
   
    # extract bam 
    def fetch(self, c=None, s=None, e=None):
        if s and e and s > e:
            (s, e) = (e, s)
        if s and s < 0:
            s = 1
        if s and s > chrLen[chrNum(c)]:
            s = chrLen[chrNum(c)] - 1
            e = chrLen[chrNum(c)] - 1
        if e and e < 0:
            s = 1
            e = 1
        if e and e > chrLen[chrNum(c)]:
            e = chrLen[chrNum(c)] - 1
        if self.downsample_ratio == 1:
            if c and s and e: 
                for a in self.bamfile.fetch(c, s, e ):
                    yield a
            else:
                for a in self.bamfile.fetch():
                    yield a
        else:
            if c and s and e:                     
                for a in  self.bamfile.fetch(c, s, e):
                    random.seed(a.query_name)
                    if random.uniform(0, 1) < self.downsample_ratio:
                        yield a
            else:                
                for a in  self.bamfile.fetch():
                    random.seed(a.query_name)
                    if random.uniform(0, 1) < self.downsample_ratio:
                        yield a
                
    # bamfile.count_coverage(contig, start, stop)

    def interval_coverage(self, i, clip=False): 
        call_args = (i[0], i[1], i[2], clip)  # i = [chrom,start,end]
        if call_args in self.interval_coverage_calls:
            return self.interval_coverage_calls[call_args]
        
        s2 = i[1]
        e2 = i[2] 
        if s2 < 0:
            s2 = 0
        if e2 > chrLen[chrNum(i[0])]:  
            e2 = chrLen[chrNum(i[0])]
        if s2 >= e2:
            return 0
        
        if clip == True or (clip is None and e2-s2 <= 1000):
            icc = sum([sum(a) for a in self.bamfile.count_coverage(i[0], s2, e2, quality_threshold=  self.mapping_quality_cutoff)]) * self.downsample_ratio / max(1.0, float(e2-s2 + 1))
            self.interval_coverage_calls[call_args] = icc
            return self.interval_coverage_calls[call_args]
        else:
            alist_len = len([a for a in self.fetch(i[0], s2, e2)
                if not a.is_unmapped and a.reference_end - 1 <= e2 and a.mapping_quality > self.mapping_quality_cutoff])
            self.interval_coverage_calls[call_args] = alist_len * self.read_length / max(1.0, float(e2 - s2 + 1))  # e2 - s2 + 1 = ws + 1
            return self.interval_coverage_calls[call_args]

    def median_coverage(self, window_size=-1, refi=-1, window_list=None): 
        if (window_size == 10000 or window_size == -1) and self.basic_stats_set and refi == -1:
            return self.downsample_stats

        num_iter = 1000
        iteri = 0
        chroffset = 0
        sumchrLen = sum([l for l in chrLen.values()]) 

        if not self.basic_stats_set:
            read_length = []
            insert_size = []
            window_list_index = 0
            non_mapping = 0
            random.seed(self.seed)

            while (window_list is not None and window_list_index < len(window_list)) or (window_list is None and iteri <= num_iter): 
                if window_list is None:
                    newpos = int(random.random() * sumchrLen) + chroffset 
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = absPos(cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000)
                if chrPos(newpos) is None:
                    logging.debug("Unable to locate reference position: ")
                    iteri+=1
                    continue

                (c, p) = chrPos(newpos)

                region = GRange([[c, p, p + 10000]], 'list')
                if c not in self.bamfile.references or p < 10000 or chrLen[chrNum(c)] < p + 10000 or len(region.intersect(conserved_regions, a_extend = 10000, b_extend = 0).gr) > 0 or len(region.intersect(centromere_list, a_extend = 10000, b_extend = 0).gr) > 0:
                    continue
                read_length += [a.infer_query_length(always=False) for a in self.fetch(c, p, p+10000) if not a.is_unmapped] 
                insert_size += [a.template_length for a in self.fetch(c, p, p+10000) if a.is_proper_pair and not a.is_reverse and a.template_length < 10000 and a.template_length > 0]
                iteri += 1

            self.read_length = np.average(read_length) 
            self.insert_size = np.average(insert_size) 
            try:
                percent_proper = len(insert_size) * 2.0 / (len(read_length) + non_mapping)
                self.percent_proper = percent_proper
            except:
                percent_proper = 1
                self.percent_proper = 1
            self.insert_std = np.std(insert_size)
            self.max_insert = self.insert_size + self.num_sdevs*self.insert_std  
            self.min_insert = max(0, self.insert_size - self.num_sdevs*self.insert_std)

        if window_size not in [-1, 300, 10000]:
            ws_list = [window_size]
        else:
            ws_list = [10000]

        wc_median = []
        wc_avg = []
        wc_std = []
        random.seed(self.seed)

        for ws in ws_list:
            wc_ws = []
            iteri = 0
            window_list_index = 0
            while (window_list is not None and window_list_index < len(window_list)) or (window_list is None and iteri <= num_iter):
                if window_list is None:
                    newpos = int(random.random() * sumchrLen) + chroffset
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = absPos( cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000 )

                if chrPos(newpos) is None:
                    logging.warning("Unable to locate reference position:")
                    iteri += 1
                    continue
                
                (c, p) = chrPos(newpos)
                region = GRange([[c, p, p + ws]], 'list')
                if c not in self.bamfile.references or p < ws or chrLen[chrNum(c)] < p + ws or len(region.intersect(conserved_regions, a_extend = 10000, b_extend = 0).gr) > 0 or len(region.intersect(centromere_list, a_extend = 10000, b_extend = 0).gr) > 0:
                    continue
                wc_ws.append(self.interval_coverage([c, p, p + ws]))
                iteri += 1
            wc_ws.sort()
            wc_ws_median = np.median(wc_ws)  
            wc_ws_filter = [c for c in wc_ws if c < 5 * wc_ws_median and c > 0]
            if len(wc_ws_filter) == 0:
                # print(len(wc_ws_filter), len(wc_ws), len([c for c in wc_ws if c > 0]), wc_ws_median)
                wc_median.append(0)
                wc_avg.append(0)
                wc_std.append(0)
            else:
                wc_median.append(wc_ws_filter[len(wc_ws_filter) // 2])
                wc_avg.append(np.average(wc_ws_filter))
                wc_std.append(np.std(wc_ws_filter))

        
        (wc_10000_median, wc_10000_avg, wc_10000_std) = (wc_median[0], wc_avg[0], wc_std[0])

        # self.pair_support = max(int(round((wc_300_avg / 10.0) * ((self.insert_size - self.read_length) / 2 / self.read_length)*self.percent_proper)), self.pair_support_min) 
        rstats = (wc_10000_median, wc_10000_avg, wc_10000_std, self.read_length,
                    self.insert_size, self.insert_std, self.min_insert, self.max_insert, 
                    self.percent_proper, self.num_sdevs, self.bamfile_filesize)
        if refi == -1:
            self.basic_stats = rstats
            self.basic_stats_set = True
            print("read length:", self.read_length, "insert size:", self.insert_size, "insert std dev:", self.insert_std,
                    "max_insert:", self.max_insert, "percent proper:", percent_proper, "num_sdevs", self.num_sdevs)
            print("coverage stats", self.basic_stats, len(wc_ws_filter))

        r = rstats
        if self.downsample < 0 or self.downsample > self.basic_stats[0]:
            self.downsample_ratio = 1
        elif self.downsample == 0:
            self.downsample_ratio = 10.0 / self.basic_stats[0] if self.basic_stats[0] > 10 else 1
        else:
            self.downsample_ratio = float(self.downsample) / self.basic_stats[0] if self.basic_stats[0] > float(self.downsample) else 1
        if self.downsample_ratio != 1:
            rr = self.downsample_ratio
            rsq = math.sqrt(rr)         
            r = [i[0] * i[1] for i in zip([rr, rr, rsq, 1, 1, 1, 1, 1, 1, 1, 1], r)]

            self.downsample_stats = r

        else:
            self.downsample_stats = self.basic_stats
