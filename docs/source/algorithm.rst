Algorithm overview
==================

ASK reconstructs amplified amplicon structures from copy-number and
breakpoint-level evidence.

1. Alignment evidence extraction
--------------------------------

ASK scans an indexed BAM file and extracts soft-clipped reads,
supplementary alignments, split-read evidence, and genomic bin counts. Reads
are filtered by mapping quality, mismatch count, duplicate status, and
alignment flags.

2. Copy-number estimation
-------------------------

Reads are counted in genomic bins. Counts are normalized and converted into
copy-number estimates. ``standard`` mode uses whole-bin counts, while ``bias``
mode uses more robust summaries for assays with uneven coverage.

3. Amplified segment detection
------------------------------

ASK identifies amplified genomic intervals from segmented copy-number
profiles. When ``--subseg`` is enabled, amplified regions can be further
refined into sub-segments.

4. Breakpoint-pair detection
----------------------------

Candidate breakpoints are inferred from clipped-read clusters and
supplementary alignments. Breakpoint pairs are filtered by read support,
distance, strand, and local copy-number context.

5. Amplicon reconstruction
--------------------------

ASK builds a graph from amplified segments and breakpoint pairs. Circular and
linear paths are reconstructed from the graph and annotated with genes, cancer
genes, and super-enhancers.

6. Targeted search
------------------

``ask-search`` follows the same evidence model but starts from a known ecDNA
structure. It limits the initial search to relevant chromosomes and breakpoint
neighborhoods, then reports matched evidence and JCS.
