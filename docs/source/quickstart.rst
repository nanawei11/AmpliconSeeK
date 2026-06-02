Quick Start
===========

Run the Example BAM
-------------------

The repository includes a small example BAM file under ``exampledata/``.
From the repository root, run:

.. code-block:: bash

   ask \
     -i exampledata/testdata.bam \
     -o exampledata/testdata \
     -g hg38 \
     --knn 3 \
     --subseg \
     --method both \
     --juncread 5 \
     --SA_with_nm \
     --run_from_pdat 1

The main circular amplicon plot is expected at:

.. code-block:: text

   exampledata/testdata_ask_plot/circular_circ_6.pdf

Run from Your Own BAM
---------------------

For a sorted and indexed BAM file:

.. code-block:: bash

   ask \
     -i sample.bam \
     -o results/sample \
     -g hg38

Prepare a BAM from FASTQ
------------------------

Paired-end reads:

.. code-block:: bash

   bwa mem -t 5 <bwa_index> test_R1.fastq.gz test_R2.fastq.gz \
     | samtools view -Shb - > test_unsorted.bam

Single-end reads:

.. code-block:: bash

   bwa mem -t 5 <bwa_index> test.fastq.gz \
     | samtools view -Shb - > test_unsorted.bam

Sort, mark duplicates, and index:

.. code-block:: bash

   samtools fixmate --threads 5 -m test_unsorted.bam - \
     | samtools sort --threads 5 -T ./ - \
     | samtools markdup --threads 5 -T ./ -S -s - test.bam

   samtools index test.bam

Common Checks
-------------

Before running ASK, confirm:

* The BAM file is coordinate-sorted.
* ``sample.bam.bai`` exists.
* The genome build matches ``-g hg19`` or ``-g hg38``.
* Soft-clipped reads and SA tags were not removed during preprocessing.
