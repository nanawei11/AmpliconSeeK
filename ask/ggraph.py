################################################################################
# GGraph class: search for circular amplicons
################################################################################
import gzip
from collections import Counter, defaultdict
from itertools import compress
import numpy as np
import pandas as pd
from memory_profiler import profile

import warnings

warnings.filterwarnings('ignore')
import misc

#------------------------------------------------------------------------------#
# GGraph class: search for circular amplicons
#------------------------------------------------------------------------------#
class GGraph:
    """
    store and analyze genomic intervels
    """

    # all nodes must be 0, 1, 2, 3, 4,...
    #--------------------------------------------------------------------------#
    def __init__(self):
        self.graph = defaultdict(list)
        self.vertex = defaultdict(list)
        self.edges = defaultdict(list)
        self.node_dict = {}
        self.single_segment_loop = []

    #--------------------------------------------------------------------------#
    def add_edge(self, a, b, weight):
        if b not in self.graph[a]:
            self.graph[a].append(b)
        if a not in self.graph[b]:
            self.graph[b].append(a)
        self.edges[(a, b)] = weight
        self.edges[(b, a)] = weight

    #--------------------------------------------------------------------------#
    def add_vertex(self, vertex_key, vertex_value):
        self.vertex[vertex_key] = vertex_value

    #--------------------------------------------------------------------------#
    def get_nodes(self):
        return list(self.graph.keys())

    #--------------------------------------------------------------------------#
    def dfs_search_circle(self, v, visited, parent, path, allpath, nodebefore):
        """
        depth-first search for circles in a graph

        Note: make visited and path have the same nodes
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        path.append(v)

        # loop for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            if v_adj not in nodebefore:
                # edge type is different from the last one
                # one bp_pair edge and one segment edge
                if (parent == -1 \
                    or self.edges[(parent, v)][0] != self.edges[(v, v_adj)][0]):
                    if (visited[v_adj] == False): # recur if the node is not visited
                        self.dfs_search_circle(v_adj, visited, v, path, allpath, nodebefore)
                    # cycle detected if the adjacent node is visited
                    # and is not the parent of current node
                    elif (parent != -1 and parent != v_adj): #and v_adj in path
                        path.append(v_adj) # add v_adj, it's the loop start
                        # store when loop detected
                        allpath.append(tuple(path)) # must use tuple! or can't run
                        path.pop() # remove v_adj from path
                        # file3 = '/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_1020/ask3_both_knn3_mr10_GBM39_1020_nofilter_circ/all_path.txt'
                        # with open(file3,'a') as f:
                        #     f.write(f'circ\t{tuple(path)}\n')
                else:
                    if (len(self.graph[v]) == 1): # end of a path, linear amplicon
                        # store when to the branch end
                        allpath.append(tuple(path)) # must use tuple! or can't run
                        # file3 = '/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_1020/ask3_both_knn3_mr10_GBM39_1020_nofilter_circ/all_path.txt'
                        # with open(file3,'a') as f:
                        #     f.write(f'linear\t{tuple(path)}\n')

        # remove current edge from path and mark it as unvisited
        path.pop()
        visited[v] = False

        return allpath

    # @profile
    #--------------------------------------------------------------------------#
    def get_all_path_contain_circle(self):
        """
        get all paths containing circle
        also include paths not containing any circle
        """

        # init onjects
        path = []
        allpath = []
        node = sorted(self.get_nodes())

        # mark all nodes as not vsisited
        visited = [False]*(max(node) + 1)

        # record nodes has already detected circle
        nodebefore = []
        # detect circle in different subgraph
        # say until all nodes are visisted
        for i in node:
            if visited[i] == False:
                self.dfs_search_circle(i, visited, -1, path, allpath, nodebefore)
            nodebefore.append(i)
        return allpath


    #--------------------------------------------------------------------------#
    # amplicon processing functions
    #--------------------------------------------------------------------------#
    def path_unique(self, circ, type = 'circular'):
        """
        get unique paths of given list of paths

        type: circular or linear
        """
        # get all circle and unique
        if (type == 'circular'):
            circ = misc.unique([row[row.index(row[-1]):-1] for row in circ])
            circ = [x for x in circ if x] # remove empty path
        elif (type == 'linear'):
            circ = [row for row in circ if (row[-1] not in row[:-1])]

        # keep ones both start and end with inner
        circ = [row for row in circ \
            if (self.edges[row[:2]][0] == 'inner' and \
                self.edges[row[-2:]][0] == 'inner')]

        if (type == 'circular'):
            # rotate path to smallest in the first
            circ = [self.rotate_path(i) for i in circ]

            # invert path to make second one smaller than last one
            circ = misc.unique([
                self.invert_path(i) if (i[1] > i[-1]) else i for i in circ])
        elif (type == 'linear'):
            circ = misc.unique([
                i[::-1] if (i[0] > i[-1]) else i for i in circ])

        return circ

    #--------------------------------------------------------------------------#
    def path_longest(self, circ, type = 'circular', check_chrom = False):
        """
        get the longest path

        if multiple longest path,
        choose the one with max edge weight

        if multiple max edge weight,
        choose the first one
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
        elif (type == 'linear'):
            last_idx = 0

        if type == 'linear':
            # get all longest
            cc = [len(row) for row in circ]
            idx = misc.all_max_index(cc)
            longest = [circ[i] for i in idx]

            # get the first one with max edge weight
            oop = []
            for row in longest:
                op = []
                for v in zip(row, row[1:] + row[:last_idx]):
                    if (self.edges[v][0] == 'outer'):
                        op.append(self.edges[v][1])
                oop.append(np.sum(op))
                circ_ref = longest[np.argmax(oop)]
        else:
            # get all longest
            cc = [len(row) for row in circ]
            circ_ref = []
            for i in set(cc):
                paths = [circ[ind] for ind, v in enumerate(cc) if v == i]
                # get the first one with max edge weight
                oop = []
                cand_ind = []
                for ind, row in enumerate(paths):
                    if check_chrom:
                        # find the all path is from same chrom 
                        chr = set([self.vertex[v][0] for v in row])
                        if len(chr) == 1:
                            cand_ind.append(ind)
                    op = []
                    for v in zip(row, row[1:] + row[:last_idx]):
                        if (self.edges[v][0] == 'outer'):
                            op.append(self.edges[v][1]) 
                    oop.append(np.sum(op))
                cand_ind.append(np.argmax(oop))
                circ_ref += [paths[t] for t in set(cand_ind)]
            # print(circ_ref)
        
        return circ_ref

    #--------------------------------------------------------------------------#
    def get_representitive_path(self, circ, type = 'circular', first_filter = False, second_filter = False, third_filter = False):
        """
        get representitive non-overlapping paths
        """
        paths = circ.copy()
        circ_repr = []
        if type == 'circular':
            t = 0
            while(paths):
                if t == 0:
                    check_chrom = True
                else:
                    check_chrom = False
                # search for longest path
                longest = self.path_longest(paths, type = type, check_chrom = check_chrom)
                circ_repr += longest
                # search non-intersected path
                # tf = [not bool(set(longest).intersection(row)) for row in paths]
                tf = []
                for row in paths:
                    tf_ = True
                    for c in longest:
                        if len(set(c).intersection(row))/len(row) > 0.7:
                            tf_ = False
                            break
                    tf.append(tf_)
                print(len(tf), sum(tf))
                # tf = [len(set(longest).intersection(row))/len(row) >= 0.4 for row in paths]
                paths = [i for (i, v) in zip(paths, tf) if v]
                t += 1
                if sum(tf) == len(paths):
                    break

        else:
            while(paths):
                # search for longest path
                longest = self.path_longest(paths, type = type)
                circ_repr.append(longest)
                # search non-intersected path
                tf = [not bool(set(longest).intersection(row)) for row in paths]
                paths = [i for (i, v) in zip(paths, tf) if v]
    
        return circ_repr

    #--------------------------------------------------------------------------#
    def make_amplicon_df(self, circ, type = 'circular'):
        """
        convert ggraph to interpretable amplicon dataframe
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
            tag = 'circ_'
        elif (type == 'linear'):
            last_idx = 0
            tag = 'line_'

        op = []
        for idx in range(len(circ)):
            row = circ[idx]
            op_seg = []
            op_bpp = []
            for v in zip(row, row[1:] + row[:last_idx]):
                if (self.edges[v][0] == 'inner'):
                    op_seg.append(self.vertex[v[0]][0:3] \
                        + self.vertex[v[1]][0:3] + self.edges[v])
                else:
                    op_bpp.append(self.edges[v])

            # put zero count for the last one of the linear amplicon
            if (type == 'linear'):
                op_bpp.append((('outer', 0)))

            for i in range(len(op_seg)):
                row_seg = op_seg[i]
                row_bpp = op_bpp[i]

                if (row_seg[2] == 'L'):
                    op.append([row_seg[0], row_seg[1], row_seg[4], \
                        '+', row_bpp[1], row_seg[7], tag + str(idx)])
                else:
                    op.append([row_seg[0], row_seg[4], row_seg[1], \
                        '-', row_bpp[1], row_seg[7], tag + str(idx)])

        colnames = ['Chrom', 'Start', 'End', 'Strand', \
            'SplitCount', 'CN', 'AmpliconID']
        df = pd.DataFrame(op, columns = colnames)

        # also add single segment loops
        if (type == 'circular'):
            node_in_path = [v for path in circ for v in path]
            ssl = [row for row in self.single_segment_loop \
                if (self.node_dict[(row[0], row[1], 'L')] not in node_in_path)]

            # add index
            ssl_df = []
            idx = len(circ)
            for row in ssl:
                ssl_df.append(row + ['circ_' + str(idx)])
                idx += 1

            # make dataframe
            ssl_df = pd.DataFrame(ssl_df, columns = colnames)
            df = pd.concat([df, ssl_df])

        return df
    
    # @profile
    #--------------------------------------------------------------------------#
    def build_ggraph_from_bp(self, bp_pair, bp_fine, seg):
        """
        build ggraph from breakpoint data
        """

        for row in bp_fine.itertuples():
            self.add_vertex(row[0], row[1:])

        for row in bp_fine.itertuples():
            self.node_dict[row[1:4]] = row[0]

        for row in bp_pair.itertuples():
            if row[1:4] in self.node_dict and row[4:7] in self.node_dict:
                a = self.node_dict[row[1:4]]
                b = self.node_dict[row[4:7]]
                w = ('outer', row[7])
                self.add_edge(a, b, w)

        for row in seg.itertuples():
            if (row[1], row[2], 'L') in self.node_dict and \
                (row[1], row[3], 'R') in self.node_dict:
                a = self.node_dict[(row[1], row[2], 'L')]
                b = self.node_dict[(row[1], row[3], 'R')]
                w = ('inner', row[4])
                # if already exist, it's a single segment loop
                if (a, b) in self.edges:
                    self.single_segment_loop.append([row[1], row[2], row[3], \
                        '+', self.edges[(a, b)][1], row[4]])
                    # don't add inner in this case
                else:
                    self.add_edge(a, b, w)


    #--------------------------------------------------------------------------#
    @staticmethod
    def invert_path(path):
        path = path[::-1]
        return path[-1:] + path[:-1]

    @staticmethod
    def rotate_path(path):
        i = path.index(min(path))
        return path[i:]+path[:i]

    @staticmethod
    def is_new_path(path, paths):
        return not path in paths


    #--------------------------------------------------------------------------#
    # currently not in use
    #--------------------------------------------------------------------------#
    def dfs_connected(self, v, visited, visited_list):
        """
        depth-first search for connected nodes (not in use currently)
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        visited_list.append(v)

        # recur for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            if (visited[v_adj] == False): # add to list if not visited
                visited_list = self.dfs_connected(v_adj, visited, visited_list)

        return visited_list

    #--------------------------------------------------------------------------#
    def get_connected_subgraph(self):
        """
        get connected subgraph nodes (not in use currently)
        """

        # init objects
        node = sorted(self.get_nodes())
        visited = []
        op = []

        # mark all nodes as not vsisited
        visited = [False]*(max(node) + 1)

        # loop to search for all subgraphs
        for v in node:
            if (visited[v] == False):
                visited_list = []
                op.append(self.dfs_connected(v, visited, visited_list))
        return op



def stat_by_seg(circ_anno, bin_norm, binsize):
    """
    Calculate copy-number statistics for each circular amplicon segment.

    For each segment in each AmpliconID, this function calculates:
      - Mean and standard deviation of CN inside the segment.
      - Mean and standard deviation of CN in the nearby left region.
      - Mean and standard deviation of CN in the nearby right region.

    Nearby left/right regions are defined by the intervals between adjacent
    segments when available. For the first or last segment on a chromosome,
    flanking windows based on `binsize` are used instead. Repeated genomic
    intervals are cached in `seg_stats` to avoid recalculating CN summaries.

    Parameters
    ----------
    circ_anno : pandas.DataFrame
        Segment-level circular amplicon annotation table. Expected columns
        include AmpliconID, Chrom, Start, End, CN, and SplitCount.
    bin_norm : pandas.DataFrame
        Bin-level normalized CN table with Chrom, Coord, and CN columns.
    binsize : int
        Bin size used to define upstream/downstream nearby regions.

    Returns
    -------
    pandas.DataFrame
        `circ_anno` with added segment CN statistics, left/right fold-change
        features, and inverse CN CV (`invCV`). Returns an empty DataFrame if
        segment statistics cannot be calculated.
    """
    # result_df = pd.DataFrame(columns=['AmpliconID', 'FragmentStart', 'FragmentEnd', 'MeanLeft', 'StdLeft', 'MeanRight', 'StdRight', 'MeanBetween', 'StdBetween'])
    seg_stats = {}
    result_dict_1 = {}
    # for each circ
    AmpliconIDs = circ_anno['AmpliconID'].unique()
    for amplicon_id in AmpliconIDs:

        # get the all segment for current circ
        fragments = circ_anno[circ_anno['AmpliconID'] == amplicon_id]
        for i in range(len(fragments)):
            #  region of current segment
            chrom = fragments.iloc[i]['Chrom']
            start = fragments.iloc[i]['Start']
            end = fragments.iloc[i]['End']
            if (chrom, start, end) in seg_stats:
                exact_mean, exact_std = seg_stats[(chrom, start, end)]
            else:
                if end - start <= binsize:
                    exact_mean = fragments.iloc[i]['CN']
                    exact_std = -1
                else:
                    exact_region = bin_norm[(bin_norm['Chrom'] == chrom) & (bin_norm['Coord'] >= start) & (bin_norm['Coord'] <= end)]
                    exact_mean = exact_region['CN'].mean()
                    if exact_region.shape[0] == 1:
                        exact_std = -1
                    else:
                        exact_std = exact_region['CN'].std()
            seg_stats[(chrom, start, end)] = (exact_mean, exact_std)   

            # left region of current segment 
            if (i == 0) or (fragments.iloc[i]['Chrom'] != fragments.iloc[i - 1]['Chrom']):
                start = fragments.iloc[i]['Start'] - 10 * binsize
                end = fragments.iloc[i]['Start'] - binsize
            else:
                start = fragments.iloc[i-1]['End']
                end = fragments.iloc[i]['Start']
                if end - start <= binsize:
                    start = start - binsize
                    end = end + binsize

            ## for ovelap seg or reback seg
            if start > end:
                start,end = end, start
                
            if (chrom, start, end) in seg_stats:
                mean_left, std_left = seg_stats[(chrom, start, end)]
            else:
                between_region = bin_norm[(bin_norm['Chrom'] == chrom) & (bin_norm['Coord'] > start) & (bin_norm['Coord'] < end)]
                mean_left = between_region['CN'].mean()
                if between_region.shape[0] == 1:
                    std_left = -1
                else:
                    std_left = between_region['CN'].std()

            seg_stats[(chrom, start, end)] = (mean_left, std_left)


            # right region of current segment 
            if (i == len(fragments)-1) or (fragments.iloc[i]['Chrom'] != fragments.iloc[i + 1]['Chrom']):
                start = fragments.iloc[i]['End'] + binsize
                end = fragments.iloc[i]['End'] + 10 * binsize
            else:
                start = fragments.iloc[i]['End']
                end = fragments.iloc[i + 1]['Start']
                if end - start <= binsize:
                    start = start - binsize
                    end = end + binsize
            
            ## for ovelap seg or reback seg
            if start > end:
                start,end = end, start
                    
            if (chrom, start, end) in seg_stats:
                mean_right, std_right = seg_stats[(chrom, start, end)]
            else:
                between_region = bin_norm[(bin_norm['Chrom'] == chrom) & (bin_norm['Coord'] > start) & (bin_norm['Coord'] < end)]
                mean_right = between_region['CN'].mean()
                if between_region.shape[0] == 1:
                    std_right = -1
                else:
                    std_right = between_region['CN'].std()

            seg_stats[(chrom, start, end)] = (mean_right, std_right)

            result_dict_1[tuple(fragments.iloc[i][:7])] = (exact_mean, exact_std, mean_left, std_left, mean_right, std_right)
    result_bin_norm = pd.DataFrame(result_dict_1).T
    if not result_bin_norm.empty:
        result_bin_norm.columns = ['Mean_Seg', 'Std_Seg', 'Left_Mean_Nearby', 'Left_Std_Nearby','Right_Mean_Nearby', 'Right_Std_Nearby']
        result_stat_df = pd.concat([circ_anno.reset_index(drop=True),result_bin_norm.reset_index(drop=True)],axis=1)
        
        result_stat_df['Mean_Seg'][result_stat_df['Mean_Seg']<0]  = 2 # advoid chrM to -999999
        result_stat_df['Mean_Seg'][result_stat_df['Mean_Seg']==0] +=  0.01
        result_stat_df['Left_Mean_Nearby'][result_stat_df['Left_Mean_Nearby']==0] +=  0.01
        result_stat_df['Right_Mean_Nearby'][result_stat_df['Right_Mean_Nearby']==0] +=  0.01
        result_stat_df['Std_Seg'][result_stat_df['Std_Seg'] == 0] +=  0.01 
        
        result_stat_df['FoldChange_left'] = np.abs(np.log2(result_stat_df['Mean_Seg']/result_stat_df['Left_Mean_Nearby']))
        result_stat_df['FoldChange_right'] =  np.abs(np.log2(result_stat_df['Mean_Seg']/result_stat_df['Right_Mean_Nearby']))
        result_stat_df['invCV'] = result_stat_df['Mean_Seg']/result_stat_df['Std_Seg']
    else:
        result_stat_df = pd.DataFrame()
    return result_stat_df


def add_stats_circ(circ_anno, bin_norm, binsize=10000):
    """
    Summarize segment-level circular amplicons and calculate the final Score.

    This function first calls `stat_by_seg` to add CN and fold-change features
    for each segment, then groups segments by AmpliconID to generate circ-level
    summary statistics. The final `Score` is based on:
      - FCleft_mean_1 and FCright_mean_1: mean left/right fold-change per segment.
      - invCNCV_mean: mean inverse CN coefficient-of-variation feature.
      - invSplitCV: inverse SplitCount coefficient-of-variation feature.

    Three penalties are applied to the base score:
      - Region imbalance penalty: absolute difference between left/right mean FC.
      - CNCV penalty: applied when invCNCV_mean is below 0.5.
      - CN_std penalty: smooth sigmoid penalty for high CN standard deviation.

    Parameters
    ----------
    circ_anno : pandas.DataFrame
        Segment-level circular amplicon annotation table.
    bin_norm : pandas.DataFrame
        Bin-level normalized CN table used by `stat_by_seg`.
    binsize : int, default 10000
        Bin size used to define nearby CN regions.

    Returns
    -------
    pandas.DataFrame
        One row per AmpliconID with circ-level summary statistics and final
        `Score`.
    """

    circ_anno = stat_by_seg(circ_anno, bin_norm, binsize)

    # the length of each segs  
    circ_anno['Length'] = circ_anno['End'] - circ_anno['Start']
    circ_anno.loc[circ_anno['invCV'] < 0, 'invCV'] = abs(circ_anno.loc[circ_anno['invCV'] < 0, 'invCV']) / 100
    # circ_anno['invCV']  = abs(circ_anno['invCV'])
    circ_anno['invCV'] = np.log(circ_anno['invCV'] + 1)
    # group by AmpliconID
    grouped = circ_anno.groupby('AmpliconID')
    
    # add stats
    result = grouped.agg({
        'Chrom': 'count', 
        'Length': 'sum',
        'SplitCount': ['sum','mean','std'],
        'CN': ['sum','mean','std'], 
        'FoldChange_left': ['sum'],
        'FoldChange_right': ['sum'],
        'invCV': ['sum','mean']
    }).reset_index()

    if 'Gene' in circ_anno.columns:
        try:
            result['Gene_num'] = grouped['Gene'].apply(lambda x: x.str.split(';').str.len().sum()).tolist()
        except:
            result['Gene_num']=0
    else:
        result['Gene_num'] = 0

    if 'CancerGene' in circ_anno.columns:  
        try:
            result['CancerGene'] = grouped['CancerGene'].apply(lambda x: x.str.split(';').str.len().sum()).tolist()
        except:
            result['CancerGene'] = 0
    else:
        result['CancerGene'] = 0
        
    if 'SE' in circ_anno.columns:
        try:
            result['SE_count'] = grouped['SE'].apply(lambda x: x.str.split(';').str.len().sum()).tolist()
        except:
            result['SE_count'] = 0
    else:
        result['SE_count'] = 0
        
    # get start and end position
    first_seg = grouped.first()[['Chrom','Start']].reset_index(drop=True)
    last_seg = grouped.last()[['Chrom','End']].reset_index(drop=True)
    result = pd.concat([result, first_seg, last_seg], axis = 1)

    result.columns = ['AmpliconID', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 'SplitCount_std', \
                      'CN_sum', 'CN_mean', 'CN_std', 'FCleft_sum', 'FCright_sum', 'invCNCV_sum', 'invCNCV_mean', 'Gene_num', \
                        'Cancergene_num', 'SE_num', 'Chrom1','Start', 'Chrom2', 'End']
    # Use -1 as a temporary sentinel for missing values before fixing std columns.
    result = result.fillna(-1)
    result.loc[result['SplitCount_std']==0,'SplitCount_std'] = 0.2
    result['invSplitCV'] = result['SplitCount_mean'] / result['SplitCount_std']
    result.loc[result['invSplitCV'] < 0, 'invSplitCV'] = abs(result.loc[result['invSplitCV'] < 0, 'invSplitCV']) / 100
    result['invSplitCV'] = np.log(result['invSplitCV'] + 1)
    result.loc[result['SplitCount_std'] == -1,'SplitCount_std'] = 1
    result.loc[result['CN_std'] == -1,'CN_std'] = 1

    # Use fold-change means derived from sum/Seg_num for a consistent score definition.
    result['FCright_mean_1'] = result['FCright_sum']/result['Seg_num']
    result['FCleft_mean_1'] = result['FCleft_sum']/result['Seg_num']
    base_score = result['FCleft_mean_1'] + result['FCright_mean_1'] + result['invCNCV_mean'] + result['invSplitCV']
    penalty_region_add = abs(result['FCleft_mean_1'] - result['FCright_mean_1'])
    penalty_CNCV = np.where(result['invCNCV_mean'] < 0.5, 1 - result['invCNCV_mean'], 0)
    penalty_CNstd = (
        0.6
        / (1 + np.exp(-0.02 * (result['CN_std'] - 200)))
        * (1 / (1 + np.exp(-0.1 * (result['CN_std'] - 50))))
    )
    scaled_penalty = np.minimum(penalty_CNCV + penalty_CNstd, 0.6)
    result['Score'] = (base_score - penalty_region_add) * (1 - scaled_penalty)
    rearrange = ['AmpliconID', 'Chrom1','Start', 'Chrom2', 'End', 'Seg_num', 'Length', 'SplitCount_sum', 
                 'SplitCount_mean', 'SplitCount_std', 'CN_sum', 'CN_mean', 'CN_std', 
                 'FCleft_sum', 'FCright_sum', 'invCNCV_sum', 'invCNCV_mean', 'invSplitCV', 
                 'Gene_num', 'Cancergene_num', 'SE_num', 'FCleft_mean_1', 'FCright_mean_1', 'Score']
    return result[rearrange]
