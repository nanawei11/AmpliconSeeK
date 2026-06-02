FAQ
===

Why is the result empty?
------------------------

Common reasons include:

* The sample does not contain a strong focal amplification.
* The sequencing depth is too low.
* ``--mincn`` is too high.
* The BAM genome build does not match ``-g``.
* Soft-clipped reads or SA tags were removed during preprocessing.
* A biased assay may need ``-d bias``.

Which data files are required?
------------------------------

For the default ``hg38`` workflow, ASK uses:

.. code-block:: text

   data/hg38_blacklist.bed
   data/hg38_refgene_process.bed12
   data/se_hg38_sort.bed
   data/hg38.genome
   data/hg38_bias.bed.gz
   data/Census_all_20200624_14_22_39.tsv
   data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai
   data/GRCh38_centromere.bed
   data/conserved_gain5_hg38.bed

Do I need the large mappability bedgraph?
-----------------------------------------

No. The current code uses ``hg38_bias.bed.gz`` for bias correction. The large
``hg38full_k35_noMM.mappability.bedgraph`` file is not directly required by
the current ASK workflow.

Can I delete ``*_ask_step*.pdat`` files?
----------------------------------------

Yes, after the run is complete. Keep them if you plan to rerun later stages
with ``--run_from_pdat``.

Why does ASK-search differ from de novo ASK?
--------------------------------------------

De novo ASK scans the BAM more broadly. ASK-search constrains evidence
collection using a known ecDNA structure, then reconstructs structures
supported in the query BAM. The output should be interpreted as targeted
evidence for the known or related structure.
