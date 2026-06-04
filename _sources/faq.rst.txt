FAQ
===

Why is the result empty?
------------------------

Common reasons include:

* The sample does not contain a strong focal amplification.
* Sequencing depth is too low.
* ``--mincn`` or ``--juncread`` is too strict.
* The BAM genome build does not match ``-g``.
* Soft-clipped reads or SA tags were removed during preprocessing.
* A coverage-biased assay may need ``-d bias``.

Which segmentation mode should I use?
-------------------------------------

Use ``standard`` for WGS and low-bias whole-genome-like data. Use ``bias`` for
ATAC-seq, scATAC-seq, ChIP-seq, WES, MNase-seq, and target-capture data.

Can I delete ``*_ask_step*.pdat`` files?
----------------------------------------

Yes, after the run is complete. Keep them if you plan to rerun later stages
with ``--run_from_pdat``.

Why does ``ask-search`` differ from de novo ``ask``?
----------------------------------------------------

De novo ASK scans the BAM more broadly. ``ask-search`` starts from a
known ecDNA structure and constrains evidence collection to relevant regions.
The output should be interpreted as targeted evidence for the known or related
structure.

Why are single-cell matrix files missing?
-----------------------------------------

ASK writes single-cell matrices only when valid cell barcodes are detected in
breakpoint-supporting reads. Bulk data and BAM files without barcode tags will
not produce these matrices.

Do I need to provide gene annotations?
--------------------------------------

Bundled annotation files are used automatically for ``hg19`` and ``hg38``.
Custom files can be supplied with ``--genefile``, ``--cgfile``, and
``--sefile`` when needed.
