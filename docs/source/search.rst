Targeted ecDNA Search
=====================

``ask-search`` tests whether a known ecDNA structure is supported in a query
BAM file. It uses the known structure to seed breakpoint-pair evidence
collection, then runs ASK-style matching, junction sequence extraction, copy
number estimation, circular reconstruction, plotting, and Junction Concordance
Score (JCS) reporting.

Run Search on the Example Data
------------------------------

First run the ``ask`` quick-start command to generate:

.. code-block:: text

   exampledata/testdata_ask_amplicon_circular_new.tsv

Then use that circular amplicon table as the known structure:

.. code-block:: bash

   ask-search \
     --circular query_sample=exampledata/testdata_ask_amplicon_circular_new.tsv \
     --bam exampledata/testdata.bam \
     --genome hg38 \
     --min-junc-cnt 5 \
     -o exampledata/testdata_search \
     --outprefix exampledata/testdata_search/testdata_search

Known Structure Format
----------------------

The ``--circular`` argument uses:

.. code-block:: text

   sample_id=known_ecDNA.tsv

ASK circular tables can be used directly. A minimal manually prepared table
should contain:

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
     - 55194959
     - +
   * - circ_0
     - chr7
     - 55222714
     - 55971771
     - +

Junction Concordance Score
--------------------------

For each reference circle, ASK-search reports:

.. code-block:: text

   JCS = validated_junctions / total_reference_junctions

A junction is considered validated when enough supporting reads match the
reference junction sequence or breakpoint-pair evidence. The default detection
threshold is ``JCS >= 0.5`` with at least five supporting reads per junction.
