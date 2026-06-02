AmpliconSeeK Documentation
==========================

AmpliconSeeK (ASK) is a Python toolkit for detecting amplified genomic
structures and candidate extrachromosomal DNA (ecDNA) from indexed alignment
files. ASK supports de novo amplicon discovery and targeted search of known
ecDNA structures in new BAM files.

The main workflows are:

* ``ask``: de novo detection from a sorted and indexed BAM file.
* ``ask-search``: targeted evidence search using a known ecDNA structure.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   search
   outputs

.. toctree::
   :maxdepth: 2
   :caption: Reference

   cli_reference
   algorithm
   faq
