Output Files
============

ASK writes files using:

.. code-block:: text

   {outprefix}_ask_{result_name}

Main De Novo Outputs
--------------------

.. list-table::
   :header-rows: 1

   * - File
     - Description
   * - ``*_ask_stats.tsv``
     - Run-level summary.
   * - ``*_ask_bin_count.tsv``
     - Raw genomic bin counts.
   * - ``*_ask_bin_count_norm.tsv``
     - Normalized bin counts and copy number estimates.
   * - ``*_ask_cn_segmentation.tsv``
     - Copy number segmentation.
   * - ``*_ask_amplified_segment.tsv``
     - Amplified genomic segments.
   * - ``*_ask_breakpoint.tsv``
     - Candidate breakpoint positions.
   * - ``*_ask_breakpoint_pair_raw.tsv``
     - Raw breakpoint-pair candidates.
   * - ``*_ask_breakpoint_pair.tsv``
     - Filtered breakpoint pairs used for reconstruction.
   * - ``*_ask_amplicon_circular.tsv``
     - Candidate circular amplicon structures.
   * - ``*_ask_amplicon_circular_new.tsv``
     - Annotated and filtered circular amplicon table.
   * - ``*_ask_amplicon_circular_stat_new.tsv``
     - Summary statistics for circular amplicons.
   * - ``*_ask_amplicon_linear.tsv``
     - Candidate linear amplicon structures.
   * - ``*_ask_alignment_sequence.tsv``
     - Read-level alignment sequence evidence.
   * - ``*_ask_junctionseq/``
     - Per-amplicon junction sequence files.
   * - ``*_ask_plot/``
     - Amplicon visualization PDFs.
   * - ``*_ask_step1.pdat`` to ``*_ask_step4.pdat``
     - Intermediate cache files for rerunning later steps.

ASK-search Outputs
------------------

In addition to ASK-style outputs, ``ask-search`` writes:

.. list-table::
   :header-rows: 1

   * - File
     - Description
   * - ``known_ecDNA_segments.tsv``
     - Parsed known ecDNA segments.
   * - ``known_ecDNA_breakpoint_pairs.tsv``
     - Reference breakpoint pairs derived from the known structure.
   * - ``known_breakpoint_seed.tsv``
     - Targeted breakpoint seed regions.
   * - ``*_ask_jcs.tsv``
     - Junction Concordance Score summary.
