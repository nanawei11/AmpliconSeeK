Targeted ecDNA search
=====================

``ask-search`` tests whether a known ecDNA structure is supported in a new
BAM file. It is useful when the ecDNA structure was discovered in one
sample and you want to search for the same or related structure in another
sample.

Basic command
-------------

.. code-block:: bash

   ask-search \
     --circular query_sample=known_ecDNA.tsv \
     --bam query_sample.bam \
     --genome hg38 \
     --min-junc-cnt 5 \
     -o search_out \
     --outprefix search_out/query_sample_search

Use ASK output as the known structure
-------------------------------------

The circular amplicon table from ``ask`` can be used directly:

.. code-block:: bash

   ask-search \
     --circular query_sample=results/sample_ask_amplicon_circular.tsv \
     --bam query_sample.bam \
     --genome hg38 \
     --min-junc-cnt 5 \
     -o search_out \
     --outprefix search_out/query_sample_search

Known structure format
----------------------

The ``--circular`` argument uses:

.. code-block:: text

   sample_id=/path/to/known_ecDNA.tsv

At minimum, the known-structure table should contain:

.. list-table::
   :header-rows: 1

   * - AmpliconID
     - Chrom
     - Start
     - End
   * - circ_0
     - chr7
     - 54830975
     - 56117062

If the segment order and strand are known, include ``Strand``:

.. list-table::
   :header-rows: 1

   * - AmpliconID
     - Chrom
     - Start
     - End
     - Strand
   * - circ_0
     - chr7
     - 54830975
     - 55200000
     - +
   * - circ_0
     - chr7
     - 55500000
     - 56117062
     - +

What search mode does
---------------------

``ask-search``:

1. Parses the known ecDNA structure.
2. Derives reference breakpoint pairs.
3. Collects breakpoint evidence only from relevant chromosomes and breakpoint
   windows.
4. Matches observed breakpoint pairs to the reference structure.
5. Reconstructs supported circular structures and writes ASK-style outputs.
6. Reports Junction Concordance Score (JCS).

Junction concordance score
--------------------------

For each reference circle:

.. code-block:: text

   JCS = validated_junctions / total_reference_junctions

By default, a reference junction is validated when it has at least five
supporting reads, and a circle is detected when ``JCS > 0.5``.

Single-cell output
------------------

If cell barcodes are detected in breakpoint-supporting reads, ``ask-search``
also writes:

.. code-block:: text

   *_ask_sc_support_matrix.tsv
   *_ask_sc_normal_alignment_matrix.tsv

Rows are breakpoint pairs and columns are cell barcodes.
