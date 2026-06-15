Overview
========

ASK is designed for sequencing assays where amplified copy-number signal and
breakpoint-level evidence can be extracted from aligned reads. It can be used
with WGS, WES, ChIP-seq, MNase-seq, ATAC-seq, scATAC-seq, and target-capture
sequencing data.

What ASK does
-------------

ASK integrates three evidence layers:

* Copy-number enrichment from genomic bin counts.
* Breakpoint evidence from soft-clipped reads, split reads, and supplementary
  alignments.
* Graph-based reconstruction of circular and linear amplicon structures.

The reconstructed amplicons are annotated with genes, cancer genes, and
super-enhancers, then visualized with ASK-style amplicon plots.

Main outputs
------------

The most commonly inspected outputs are:

.. list-table::
   :header-rows: 1

   * - File
     - Meaning
   * - ``*_ask_amplicon_circular.tsv``
     - Candidate circular amplicons.
   * - ``*_ask_breakpoint_pair.tsv``
     - Breakpoint pairs used for reconstruction.
   * - ``*_ask_amplicon_circular_stat.tsv``
     - Summary and score for each circular amplicon.
   * - ``*_ask_plot/``
     - PDF visualization of amplified segments and amplicons.
   * - ``*_ask_sc_support_matrix.tsv``
     - Single-cell junction-support matrix, generated only when cell barcodes
       are detected.

Recommended reading order
-------------------------

New users should start with :doc:`installation`, then run the example command
in :doc:`quickstart`. After that, use :doc:`outputs` to interpret the result
files. Users who already have a known ecDNA structure can use the ``ask-search`` command described in
:doc:`cli_reference`.
