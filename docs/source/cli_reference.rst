Command Line Reference
======================

De Novo Detection
-----------------

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
     - Input BAM file.
   * - ``-o, --outprefix``
     - input BAM prefix
     - Output prefix.
   * - ``-g, --genome``
     - ``hg38``
     - Genome build, usually ``hg19`` or ``hg38``.
   * - ``-w, --binsize``
     - ``10000``
     - Bin size for copy number estimation.
   * - ``-k, --mapq``
     - ``20``
     - Minimum mapping quality.
   * - ``-l, --nmmax``
     - ``1``
     - Maximum NM mismatches.
   * - ``-c, --mincn``
     - ``5``
     - Minimum copy number for amplified segments.
   * - ``-d, --segmode``
     - ``standard``
     - Segmentation mode: ``standard`` or ``bias``.
   * - ``-n, --bpcount``
     - ``5``
     - Minimum clipped reads for candidate breakpoints.
   * - ``-m, --juncread``
     - ``10``
     - Minimum reads for breakpoint pairs.
   * - ``--SA_with_nm``
     - ``False``
     - Use NM mismatch values when filtering supplementary alignments.
   * - ``--subseg``
     - ``False``
     - Infer sub-segments within amplified regions.
   * - ``--method``
     - ``both``
     - Breakpoint-pair strategy.
   * - ``--knn``
     - ``3``
     - Neighbor parameter used in segment graph construction.
   * - ``--run_from_pdat``
     - ``0``
     - Resume from saved intermediate ``.pdat`` files.

Targeted Search
---------------

.. code-block:: bash

   ask-search \
     --circular query_sample=known_ecDNA.tsv \
     --bam sample.bam \
     --genome hg38 \
     -o search_out \
     --outprefix search_out/sample_search

Common options:

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
   * - ``--genome``
     - ``hg38``
     - Genome build for default annotation files.
   * - ``--window``
     - ``200``
     - Window around known breakpoints for targeted evidence collection.
   * - ``--mapq``
     - ``20``
     - Minimum mapping quality.
   * - ``--nmmax``
     - ``1``
     - Maximum NM mismatches.
   * - ``--min-junc-cnt``
     - ``5``
     - Minimum junction support count.
   * - ``--bpp-min-dist``
     - ``50``
     - Minimum distance between breakpoint-pair positions.
   * - ``--jcs-threshold``
     - ``0.5``
     - JCS threshold for detection.
   * - ``--jcs-min-support``
     - ``5``
     - Minimum supporting reads for a validated junction.
