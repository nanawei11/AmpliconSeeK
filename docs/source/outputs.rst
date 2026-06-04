Output files
============

ASK writes files using the output prefix:

.. code-block:: text

   {outprefix}_ask_{result_name}

For example, ``-o results/sample`` creates files such as
``results/sample_ask_amplicon_circular.tsv``.

Core result files
-----------------

.. list-table::
   :header-rows: 1

   * - File or directory
     - Description
   * - ``*_ask_amplicon_circular.tsv``
     - Candidate circular amplicon structures.
   * - ``*_ask_amplicon_circular_stat.tsv``
     - Summary statistics and score for each circular amplicon.
   * - ``*_ask_amplicon_linear.tsv``
     - Candidate linear amplicon structures.
   * - ``*_ask_amplified_segment.tsv``
     - Amplified genomic segments inferred from copy-number signal.
   * - ``*_ask_breakpoint.tsv``
     - Candidate breakpoint positions.
   * - ``*_ask_breakpoint_pair.tsv``
     - Final breakpoint pairs used for amplicon reconstruction.
   * - ``*_ask_breakpoint_pair_raw.tsv``
     - Raw breakpoint-pair candidates before final filtering.
   * - ``*_ask_bin_count.tsv``
     - Raw genomic bin counts.
   * - ``*_ask_bin_count_norm.tsv``
     - Normalized bin counts and copy-number estimates.
   * - ``*_ask_cn_segmentation.tsv``
     - Copy-number segmentation.
   * - ``*_ask_junctionseq/``
     - Per-amplicon junction sequence files, including breakpoint-pair
       junction sequences and supporting read sequences aligned across each
       junction.
   * - ``*_ask_plot/``
     - PDF figures for amplified segments and amplicons.
   * - ``*_ask_stats.tsv``
     - Plain-text run summary.

Targeted search files
---------------------

``ask-search`` additionally writes:

.. list-table::
   :header-rows: 1

   * - File
     - Description
   * - ``known_ecDNA_segments.tsv``
     - Parsed known ecDNA segments.
   * - ``known_ecDNA_breakpoint_pairs.tsv``
     - Reference breakpoint pairs derived from the known structure.
   * - ``known_breakpoint_seed.tsv``
     - Breakpoint seed regions used for targeted evidence collection.
   * - ``*_ask_jcs.tsv``
     - Junction Concordance Score summary.

Single-cell files
-----------------

For single-cell assays, ASK writes barcode-level matrices only when valid cell
barcodes are detected:

.. list-table::
   :header-rows: 1

   * - File
     - Description
   * - ``*_ask_sc_support_matrix.tsv``
     - Junction-support matrix. Rows are breakpoint pairs and columns are cell
       barcodes.
   * - ``*_ask_sc_normal_alignment_matrix.tsv``
     - Normal-alignment matrix with the same rows and columns as the support
       matrix.

Important columns
-----------------

.. list-table::
   :header-rows: 1

   * - Column
     - Where it appears
     - Meaning
   * - ``AmpliconID``
     - Amplicon tables
     - Reconstructed amplicon identifier.
   * - ``Chrom``, ``Start``, ``End``
     - Segment and amplicon tables
     - Genomic interval.
   * - ``CN``
     - Copy-number and amplicon tables
     - Copy-number estimate.
   * - ``Chrom1``, ``Coord1``, ``Clip1``
     - Breakpoint-pair tables
     - First breakpoint side.
   * - ``Chrom2``, ``Coord2``, ``Clip2``
     - Breakpoint-pair tables
     - Second breakpoint side.
   * - ``Count``
     - Breakpoint-pair tables
     - Supporting read count.
   * - ``Readbarcode``
     - Breakpoint-pair tables
     - Cell barcodes from single-cell data that support the breakpoint pair.
       Empty lists indicate that no valid cell barcode was detected for that
       breakpoint pair.
