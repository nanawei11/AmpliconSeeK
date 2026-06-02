Installation
============

Create a Conda Environment
--------------------------

Install Miniconda by following the official Conda instructions, then configure
the common bioinformatics channels:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

Create and activate an environment for ASK:

.. code-block:: bash

   conda create -n ask --no-channel-priority \
     pysam pandas numpy matplotlib statsmodels seaborn scipy scikit-learn
   conda activate ask

Install ASK from the repository root:

.. code-block:: bash

   pip install .

Check the command line entry points:

.. code-block:: bash

   ask --help
   ask-search --help

Required Input
--------------

ASK expects a sorted BAM file and an index file in the same directory:

.. code-block:: text

   sample.bam
   sample.bam.bai

The BAM should be aligned to the same genome build selected with ``-g`` or
``--genome``. The bundled annotation files currently support ``hg19`` and
``hg38``.
