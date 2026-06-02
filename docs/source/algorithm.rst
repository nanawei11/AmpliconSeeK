Algorithm Overview
==================

ASK combines copy number signal and breakpoint-level evidence to reconstruct
candidate amplified structures.

1. Alignment Processing
-----------------------

ASK scans the indexed BAM file, extracts soft-clipped reads, supplementary
alignment evidence, and genomic bin counts. Reads are filtered by mapping
quality, mismatch count, duplicate status, and alignment flags.

2. Copy Number Estimation
-------------------------

Reads are counted in genomic bins. Counts are normalized and converted into
copy number estimates. In ``bias`` mode, GC and mappability bias correction is
applied using the bundled ``*_bias.bed.gz`` annotation.

3. Amplified Segment Detection
------------------------------

Segmented copy number profiles are used to identify amplified intervals.
Segments can be further refined when ``--subseg`` is enabled.

4. Breakpoint-Pair Detection
----------------------------

Candidate breakpoints are inferred from clipped-read clusters and
supplementary alignments. Breakpoint pairs are filtered by read support,
distance, strand, and local copy number context.

5. Amplicon Reconstruction
--------------------------

ASK builds a graph from amplified segments and breakpoint pairs. Candidate
circular and linear amplicon paths are reconstructed from the graph and then
annotated with genes, cancer genes, and super-enhancers.

6. Targeted Search
------------------

ASK-search starts from a known ecDNA structure. It restricts initial evidence
collection to the chromosomes and breakpoint windows implied by the known
structure, then runs ASK-style breakpoint matching, junction sequence
extraction, reconstruction, plotting, and JCS scoring.
