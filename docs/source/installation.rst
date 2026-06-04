Installation
============

ASK can be installed from PyPI or from the GitHub source tree. A Conda
environment is recommended because ASK depends on common scientific Python and
bioinformatics packages.

Create a Conda environment
--------------------------

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

   conda create -n ask --no-channel-priority \
     pysam pandas numpy matplotlib statsmodels seaborn scipy scikit-learn
   conda activate ask

Install from PyPI
-----------------

.. code-block:: bash

   pip install ask-ecdna

Install from GitHub
-------------------

Use this option if you want the latest development version:

.. code-block:: bash

   git clone https://github.com/nanawei11/AmpliconSeeK.git
   cd AmpliconSeeK
   pip install .

Check the installation
----------------------

.. code-block:: bash

   ask --help
   ask-search --help

Required input
--------------

ASK expects a coordinate-sorted BAM file and its index:

.. code-block:: text

   sample.bam
   sample.bam.bai

The BAM genome build must match the ``-g`` or ``--genome`` option. The
bundled annotation files support ``hg19`` and ``hg38``.
