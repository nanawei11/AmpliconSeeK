Quick start
===========

This page shows the shortest path from an indexed BAM file to ASK output.

Get the example data
--------------------

The example BAM file, ``exampledata/testdata.bam``, is available in the
AmpliconSeeK GitHub repository:

https://github.com/nanawei11/AmpliconSeeK/tree/main/exampledata

If you clone the repository, the example data are already in the expected
location:

.. code-block:: bash

   git clone https://github.com/nanawei11/AmpliconSeeK.git
   cd AmpliconSeeK

Run the example BAM
-------------------

From the repository root:

.. code-block:: bash

   ask \
     -i exampledata/testdata.bam \
     -o exampledata/testdata/samplename \
     -g hg38 \
     --subseg \
     --juncread 5 \
     --SA_with_nm

The main output files are written with the prefix
``exampledata/testdata/samplename_ask``.

Run your own BAM
----------------

For WGS or other whole-genome-like data:

.. code-block:: bash

   ask \
     -i sample.bam \
     -o results/sample \
     -g hg38 \
     --subseg \
     --juncread 5 \
     --SA_with_nm

For coverage-biased assays such as ATAC-seq, scATAC-seq, ChIP-seq, WES,
MNase-seq, or target-capture sequencing, use ``-d bias``:

.. code-block:: bash

   ask \
     -i sample.bam \
     -o results/sample \
     -g hg38 \
     -d bias \
     --subseg \
     --juncread 5 \
     --SA_with_nm

Choose segmentation mode
------------------------

``-d/--segmode`` tells ASK what type of coverage profile to expect.

.. list-table::
   :header-rows: 1

   * - Mode
     - Recommended data
   * - ``standard``
     - WGS and low-bias whole-genome-like data. This is the default.
   * - ``bias``
     - ATAC-seq, scATAC-seq, ChIP-seq, WES, MNase-seq, and target-capture
       sequencing.

Prepare a BAM from FASTQ
------------------------

Paired-end FASTQ:

.. code-block:: bash

   bwa mem -t 5 <bwa_index> test_R1.fastq.gz test_R2.fastq.gz \
     | samtools view -Shb - > test_unsorted.bam

Single-end FASTQ:

.. code-block:: bash

   bwa mem -t 5 <bwa_index> test.fastq.gz \
     | samtools view -Shb - > test_unsorted.bam

Sort, mark duplicates, and index:

.. code-block:: bash

   samtools fixmate --threads 5 -m test_unsorted.bam - \
     | samtools sort --threads 5 -T ./ - \
     | samtools markdup --threads 5 -T ./ -S -s - test.bam

   samtools index test.bam

Checklist
---------

Before running ASK, confirm that:

* The BAM file is coordinate sorted.
* The BAM index exists.
* The genome build matches ``-g hg19`` or ``-g hg38``.
* Soft-clipped reads and SA tags were not removed during preprocessing.
