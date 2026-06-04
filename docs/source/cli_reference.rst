Command line reference
======================

This page lists the most commonly used options. Run ``ask --help`` or
``ask-search --help`` for the full command-line help.

``ask``
-------

Basic usage:

.. code-block:: bash

   ask -i sample.bam -o results/sample -g hg38

Core options:

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``-i, --bamfile``
     - required
     - Input sorted BAM file.
   * - ``-o, --outprefix``
     - input prefix
     - Output prefix.
   * - ``-g, --genome``
     - ``hg38``
     - Genome build, usually ``hg19`` or ``hg38``.
   * - ``-w, --binsize``
     - ``10000``
     - Bin size for copy-number estimation.
   * - ``-d, --segmode``
     - ``standard``
     - Input data mode: ``standard`` or ``bias``.
   * - ``-k, --mapq``
     - ``20``
     - Minimum mapping quality.
   * - ``-l, --nmmax``
     - ``1``
     - Maximum NM mismatch count.
   * - ``-c, --mincn``
     - ``5``
     - Minimum copy number for amplified segments.
   * - ``-n, --bpcount``
     - ``5``
     - Minimum clipped reads for candidate breakpoints.
   * - ``-m, --juncread``
     - ``10``
     - Minimum reads for breakpoint pairs.
   * - ``--subseg``
     - ``False``
     - Infer sub-segments within amplified regions.
   * - ``--SA_with_nm``
     - ``False``
     - Use NM mismatch values when filtering supplementary alignments.
   * - ``--run_from_pdat``
     - ``0``
     - Resume from saved intermediate ``.pdat`` files.

``ask-search``
--------------

Basic usage:

.. code-block:: bash

   ask-search \
     --circular query_sample=known_ecDNA.tsv \
     --bam query_sample.bam \
     -o search_out \
     --outprefix search_out/query_sample_search

Core options:

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--circular``
     - required
     - Known ecDNA structure in ``sample_id=table.tsv`` format.
   * - ``--bam``
     - required
     - Query BAM file.
   * - ``-o, --outdir``
     - required
     - Output directory.
   * - ``--outprefix``
     - ``outdir/<bam-stem>``
     - ASK-style output prefix.
   * - ``--genome``
     - ``hg38``
     - Genome build for bundled annotation files.
   * - ``--window``
     - ``200``
     - Breakpoint-neighborhood search window.
   * - ``--mapq``
     - ``20``
     - Minimum mapping quality.
   * - ``--nmmax``
     - ``1``
     - Maximum NM mismatch count.
   * - ``--min-junc-cnt``
     - ``1``
     - Minimum junction read count for circular reconstruction.
   * - ``--jcs-min-support``
     - ``5``
     - Minimum reads required to validate one reference junction.
   * - ``--min-jcs``
     - ``0.5``
     - Circle-level JCS detection threshold.
