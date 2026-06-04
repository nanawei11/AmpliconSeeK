AmpliconSeeK documentation
==========================

AmpliconSeeK (ASK) detects amplified genomic regions, breakpoint-pair
evidence, and candidate extrachromosomal DNA (ecDNA) structures from indexed
alignment files. ASK supports both de novo discovery and targeted search of
known ecDNA structures across samples.

ASK provides two command-line workflows:

.. list-table::
   :header-rows: 1

   * - Workflow
     - Command
     - Use case
   * - De novo detection
     - ``ask``
     - Start from a sorted, indexed BAM file and reconstruct candidate
       circular or linear amplicons.
   * - Targeted search
     - ``ask-search``
     - Test whether a known ecDNA structure is supported in another BAM
       file.

Quick links
-----------

* Install ASK: :doc:`installation`
* Run the example data: :doc:`quickstart`
* Search a known ecDNA structure: :doc:`search`
* Interpret output files: :doc:`outputs`
* Check command-line options: :doc:`cli_reference`

.. toctree::
   :maxdepth: 2
   :caption: User guide

   overview
   installation
   quickstart
   search
   outputs

.. toctree::
   :maxdepth: 2
   :caption: Reference

   cli_reference
   algorithm
