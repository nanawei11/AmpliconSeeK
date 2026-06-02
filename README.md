# AmpliconSeeK (ASK)

AmpliconSeeK (ASK) is a Python toolkit for detecting and reconstructing amplified genomic structures and candidate extrachromosomal DNA (ecDNA) from indexed alignment files, supporting both de novo discovery and targeted search of known ecDNA structures.

**Latest Release:**

* Github: v3

## Table of Contents

- [Overview](#overview)
- [Software dependencies](#software-dependencies)
- [Installation](#installation)
  - [How to install python and pre-required packages](#how-to-install-python-and-pre-required-packages)
- [Input Data Preparation](#input-data-preparation)
  - [Required Input Data](#required-input-data)
  - [BAM](#bam)
  - [Reference Annotation Data](#reference-annotation-data)
  - [Known ecDNA Structure for Search](#known-ecdna-structure-for-search)
- [Directory Structure](#directory-structure)
- [De Novo ecDNA Detection](#de-novo-ecdna-detection)
  - [How to run from bam file](#how-to-run-from-bam-file)
  - [Output](#output)
  - [How to prepare bam file](#how-to-prepare-bam-file)
- [Targeted ecDNA Search](#targeted-ecdna-search)
  - [What Search Mode Does](#what-search-mode-does)
  - [Parameters](#parameters)
  - [Output](#output-1)
- [Output Files](#output-files)
- [File Formats](#file-formats)
  - [Circular Amplicon File](#circular-amplicon-file)
  - [Circular Amplicon Statistics File](#circular-amplicon-statistics-file)
  - [Breakpoint-Pair File](#breakpoint-pair-file)
  - [JCS File](#jcs-file)
- [Algorithm Overview](#algorithm-overview)
- [Checkpointing and Modular Usage](#checkpointing-and-modular-usage)
- [License](#license)
- [Contact](#contact)

## Overview

Extrachromosomal DNA (ecDNA) is a dynamic form of oncogene amplification that contributes to cancer progression through high-copy gene dosage, regulatory rewiring, and cell-to-cell heterogeneity. AmpliconSeeK (ASK) is a computational framework for identifying ecDNA-associated amplicon structures from diverse high-throughput sequencing data, including WGS, WES, ChIP-seq, MNase-seq, ATAC-seq, scATAC-seq, and target-capture sequencing. ASK integrates copy-number signal from genomic bin counts with breakpoint-level evidence, including soft-clipped reads, split reads, supplementary alignments, breakpoint pairs, and junction sequences, to infer amplified segments and reconstruct candidate circular or linear amplicons. Candidate structures are annotated with genes, cancer genes, and super-enhancers and visualized with ASK-style amplicon plots.

ASK provides two main workflows:


| Workflow          | Command      | Description                                                                                             |
| ----------------- | ------------ | ------------------------------------------------------------------------------------------------------- |
| De novo detection | `ask`        | Detect amplified segments, breakpoint pairs, and candidate circular amplicons directly from a BAM file. |
| Targeted search   | `ask-search` | Search a new BAM file for evidence supporting a known ecDNA structure.                                 |

ASK can be applied to sequencing assays with genomic alignment signals, including WGS, WES, ChIP-seq, MNase-seq, ATAC-seq, scATAC-seq, and target-capture sequencing.

## Software dependencies

* The software has been tested in MacOSX and Linux system.
* The software does not depend on any other softwares except some basic python packages.
* Pre-required python packages: pysam, pandas, numpy, statsmodels, matplotlib, seaborn

## Installation

### How to install python and pre-required packages

Install Miniconda by following https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Set up bioconda channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create an environment with the pre-required python packages:

```bash
conda create -n ask --no-channel-priority pysam pandas numpy matplotlib statsmodels seaborn scipy scikit-learn
```

Activate environment and install ASK:

```basha
conda activate ask
pip install ask-ecdna
```

Now, you are ready to run ASK:

```bash
ask --help
ask-search --help
```

## Input Data Preparation

### Required Input Data

ASK requires the following input data:


| Data Type                 | Required for`ask` | Required for`ask-search` | Description                                             |
| ------------------------- | ----------------- | ------------------------ | ------------------------------------------------------- |
| BAM                       | Yes               | Yes                      | Sorted and indexed alignment file                       |
| BAM index                 | Yes               | Yes                      | `.bai`  index file                                     |
| Genome annotation         | Recommended       | Recommended              | Gene annotation BED12 file                              |
| Cancer gene list          | Optional          | Optional                 | Cancer gene census file                                 |
| Super-enhancer annotation | Optional          | Optional                 | BED file for SE annotation                              |
| Known ecDNA structure     | No                | Yes                      | ASK circular table or manually prepared known structure |

### BAM

The input alignment file should be sorted and indexed:

```text
sample.bam
sample.bam.bai
```

For de novo detection, duplicate marking is recommended before running ASK.

### Reference Annotation Data

ASK includes commonly used annotation files under `data/`. For example, with `--genome hg38`, ASK expects files such as:

```text
data/hg38_refgene_process.bed12
data/se_hg38_sort.bed
data/Census_all_20200624_14_22_39.tsv
```

Custom annotation files can be provided manually:

```bash
--genefile /path/to/gene.bed12
--sefile /path/to/super_enhancer.bed
--cgfile /path/to/cancer_gene.tsv
```

The genome build used by the BAM file and annotation files should match.

### Known ecDNA Structure for Search

`ask-search` accepts an ASK circular amplicon table:

```text
*_ask_amplicon_circular.tsv
```

It also accepts a manually prepared known-structure table. At minimum, the table should contain:


| AmpliconID | Chrom |    Start |      End |
| ---------- | ----- | -------: | -------: |
| circ_0     | chr7  | 54830975 | 56117062 |

If segment-level order and strand are available, include them:


| AmpliconID | Chrom |    Start |      End | Strand |
| ---------- | ----- | -------: | -------: | ------ |
| circ_0     | chr7  | 54830975 | 55200000 | +      |
| circ_0     | chr7  | 55500000 | 56117062 | +      |

ASK uses the known structure to derive reference breakpoint pairs for targeted search.

## Directory Structure

A typical ASK project can be organized as:

```text
/path/to/ask_project/
├── data/
│   ├── hg38_refgene_process.bed12
│   ├── se_hg38_sort.bed
│   └── Census_all_20200624_14_22_39.tsv
├── bam/
│   ├── sample.bam
│   └── sample.bam.bai
├── known_ecDNA/
│   └── known_ecDNA.tsv
├── ask_denovo/
│   ├── sample_ask_amplicon_circular.tsv
│   ├── sample_ask_breakpoint_pair.tsv
│   └── sample_ask_junctionseq/
└── ask_search/
    ├── sample_search_ask_amplicon_circular.tsv
    ├── sample_search_ask_jcs.tsv
    └── sample_search_ask_junctionseq/
```

ASK output filenames follow this convention:

```text
{outprefix}_ask_{result_name}.tsv
```

For example:

```text
sample_ask_amplicon_circular.tsv
sample_ask_breakpoint_pair.tsv
sample_ask_bin_count_norm.tsv
```

## De Novo ecDNA Detection

### How to run from bam file

Run the example BAM file included in this repository:

```bash
cd /path/to/AmpliconSeeK

ask \
  -i exampledata/testdata.bam \
  -o exampledata/testdata/samplename \
  -g hg38 \
  --subseg \
  --juncread 5 \
  --SA_with_nm
```

### Output

The command generates ASK-style output:

```text
testdata/
├──  samplename_ask_amplicon_circular.tsv
├──  samplename_ask_amplicon_circular_stat.tsv
├──  samplename_ask_amplicon_linear.tsv
├──  samplename_ask_amplified_segment.tsv
├──  samplename_ask_bin_count.tsv
├──  samplename_ask_bin_count_norm.tsv
├──  samplename_ask_breakpoint.tsv
├──  samplename_ask_breakpoint_pair.tsv
├──  samplename_ask_breakpoint_pair_raw.tsv
├──  samplename_ask_breakpoint_seg.tsv
├──  samplename_ask_clip_count.bedgraph
├──  samplename_ask_cn_segmentation.tsv
├──  samplename_ask_junctionseq
│   ├──  circ_0.tsv
│   ├──  circ_1.tsv
│   ├──  circ_2.tsv
│   └──  circ_3.tsv
├──  samplename_ask_plot
│   ├──  ampseg_0.pdf
│   ├──  circular_circ_0.pdf
│   ├──  circular_circ_1.pdf
│   ├──  circular_circ_2.pdf
│   └──  circular_circ_3.pdf
├──  samplename_ask_stats.tsv
├──  samplename_ask_step1.pdat
├──  samplename_ask_step2.pdat
├──  samplename_ask_step3.pdat
└──  samplename_ask_step4.pdat


```

### How to prepare bam file

Map FASTQ file to the genome:

```bash
# paired end
bwa_index=/path/to/hg38.fa 
bwa mem -t 5 ${bwa_index} test_R1.fastq.gz test_R2.fastq.gz | samtools view -Shb - > test_unsorted.bam

# single end
bwa mem -t 5 ${bwa_index} test.fastq.gz | samtools view -Shb - > test_unsorted.bam
```

Sort and mark duplicates:

```bash
samtools fixmate --threads 5 -m test_unsorted.bam - \
    | samtools sort --threads 5 -T ./ - \
    | samtools markdup --threads 5 -T ./ -S -s - test.bam
```

Make index:

```bash
samtools index test.bam
```

## Targeted ecDNA Search

Use a known ecDNA structure and a new BAM. For the example data, first run the `ask` command above, then use its circular amplicon table as the known structure:

```bash
ask-search \
  --circular query_sample=exampledata/testdata/samplename_ask_amplicon_circular.tsv \
  --bam exampledata/testdata.bam \
  --genome hg38 \
  --target-genes EGFR,MDM4,PDGFRA \
  --min-junc-cnt 5 \
  -o exampledata/testdata_search/testdata_search 
```

If running directly from the source tree:

```bash
python ask/ecDNA_search.py \
  --circular query_sample=exampledata/testdata/samplename_ask_amplicon_circular.tsv\
  --bam exampledata/testdata.bam \
  --genome hg38 \
  --min-junc-cnt 5 \
  -o exampledata/testdata_search/testdata_search
```

### What Search Mode Does

`ask-search` is a targeted workflow:

1. Parse the known ecDNA structure.
2. Derive reference breakpoint pairs from the known segments.
3. Collect reads around relevant chromosomes and breakpoint neighborhoods.
4. Match observed breakpoint-pair evidence to the reference breakpoint pairs.
5. Reconstruct supported circular structures from the observed evidence.
6. Report ASK-style outputs and Junction Concordance Score.

### Parameters


| Parameter           | Required | Default             | Description                                                               |
| ------------------- | -------- | ------------------- | ------------------------------------------------------------------------- |
| `--circular`        | Yes      | -                   | Known ecDNA structure in`sample_id=known_ecDNA.tsv` format                |
| `--bam`             | Yes      | -                   | Query BAM/CRAM file                                                       |
| `-o`, `--outdir`    | Yes      | -                   | Output directory                                                          |
| `--outprefix`       | No       | `outdir/<bam-stem>` | ASK-style output prefix                                                   |
| `--genome`          | No       | `hg38`              | Genome build for default annotation files                                 |
| `--target-genes`    | No       | None                | Optional comma-separated cancer genes used to filter reference structures |
| `--window`          | No       | `200`               | Breakpoint-neighborhood search window in bp                               |
| `--mapq`            | No       | `20`                | Minimum mapping quality                                                   |
| `--nmmax`           | No       | `1`                 | Maximum NM mismatch count                                                 |
| `--min-junc-cnt`    | No       | `1`                 | Minimum junction read count used before DFS circular reconstruction       |
| `--bpp-min-dist`    | No       | `50`                | Minimum same-chromosome breakpoint-pair distance in bp                    |
| `--jcs-min-support` | No       | `5`                 | Minimum supporting reads required to validate one reference junction      |
| `--min-jcs`         | No       | `0.5`               | Circle-level JCS detection threshold                                      |

### Output

The command generates ASK-style search output:

```text
ask_search/
├── known_breakpoint_seed.tsv
├── known_ecDNA_breakpoint_pairs.tsv
├── known_ecDNA_segments.tsv
├── sample_search_ask_alignment_sequence.tsv
├── sample_search_ask_amplicon_circular_new.tsv
├── sample_search_ask_amplicon_circular_stat_new.tsv
├── sample_search_ask_amplicon_linear.tsv
├── sample_search_ask_amplified_segment.tsv
├── sample_search_ask_bin_count.tsv
├── sample_search_ask_bin_count_norm.tsv
├── sample_search_ask_breakpoint_pair.tsv
├── sample_search_ask_breakpoint_pair_raw.tsv
├── sample_search_ask_breakpoint_seq.tsv
├── sample_search_ask_breakpoint.tsv
├── sample_search_ask_clip_count.bedgraph
├── sample_search_ask_cn_segmentation.tsv
├── sample_search_ask_jcs.tsv
├── sample_search_ask_stats.tsv
├── sample_search_ask_step1.pdat
├── sample_search_ask_step2.pdat
├── sample_search_ask_step3.pdat
├── sample_search_ask_step4.pdat
├── sample_search_ask_junctionseq/
└── plot/
```

## Output Files


| File or Directory                        | Generated by        | Description                                                     |
| ---------------------------------------- | ------------------- | --------------------------------------------------------------- |
| `*_ask_amplicon_circular.tsv`            | `ask`, `ask-search` | Candidate circular amplicon/ecDNA structures                    |
| `*_ask_amplicon_circular_stat.tsv`       | `ask`, `ask-search` | Summary statistics for circular amplicons                       |
| `*_ask_amplicon_linear.tsv`              | `ask`, `ask-search` | Candidate linear amplicon structures                            |
| `*_ask_amplified_segment.tsv`            | `ask`, `ask-search` | Amplified genomic segments inferred from copy number signal     |
| `*_ask_breakpoint.tsv`                   | `ask`, `ask-search` | Candidate breakpoint positions                                  |
| `*_ask_breakpoint_pair.tsv`              | `ask`, `ask-search` | Final breakpoint pairs used for amplicon reconstruction         |
| `*_ask_breakpoint_pair_raw.tsv`          | `ask`, `ask-search` | Raw breakpoint-pair candidates before final filtering           |
| `*_ask_breakpoint_seq.tsv`               | `ask`, `ask-search` | Breakpoint-associated sequence information                      |
| `*_ask_alignment_sequence.tsv`           | `ask`, `ask-search` | Read-level alignment sequence evidence for breakpoint junctions |
| `*_ask_junctionseq/`                     | `ask`, `ask-search` | Per-amplicon junction sequence files                            |
| `*_ask_bin_count.tsv`                    | `ask`, `ask-search` | Raw genomic bin counts                                          |
| `*_ask_bin_count_norm.tsv`               | `ask`, `ask-search` | Normalized bin counts for copy number estimation                |
| `*_ask_cn_segmentation.tsv`              | `ask`, `ask-search` | Copy number segmentation result                                 |
| `*_ask_clip_count.bedgraph`              | `ask`, `ask-search` | Soft-clipping evidence track                                    |
| `*_ask_stats.tsv`                        | `ask`, `ask-search` | Run-level summary statistics                                    |
| `*_ask_step1.pdat` to `*_ask_step4.pdat` | `ask`, `ask-search` | Intermediate cache files                                        |
| `*_ask_jcs.tsv`                          | `ask-search`        | Junction Concordance Score summary                              |
| `known_ecDNA_segments.tsv`               | `ask-search`        | Parsed known ecDNA segments used as the search target           |
| `known_ecDNA_breakpoint_pairs.tsv`       | `ask-search`        | Reference breakpoint pairs derived from the known structure     |
| `known_breakpoint_seed.tsv`              | `ask-search`        | Breakpoint seed table used for targeted evidence collection     |
| `plot/`                                  | `ask`, `ask-search` | Amplicon visualization figures                                  |

## File Formats

### Circular Amplicon File

`*_ask_amplicon_circular_new.tsv` reports candidate circular amplicon structures. Each row describes a genomic segment assigned to a candidate circular structure.


| Column         | Description                            |
| -------------- | -------------------------------------- |
| `AmpliconID`   | Candidate circular amplicon identifier |
| `Chrom`        | Chromosome name                        |
| `Start`, `End` | Genomic segment coordinates            |
| `Strand`       | Segment orientation when available     |
| `CN`           | Copy number estimate when available    |
| `Gene`         | Overlapping gene annotation            |
| `CancerGene`   | Overlapping cancer gene annotation     |
| `SE`           | Overlapping super-enhancer annotation  |

### Circular Amplicon Statistics File

`*_ask_amplicon_circular_stat_new.tsv` summarizes each candidate circle.


| Column                              | Description                                                    |
| ----------------------------------- | -------------------------------------------------------------- |
| `AmpliconID`                        | Candidate circular amplicon identifier                         |
| `Seg_num`                           | Number of segments in the circle                               |
| `Length`                            | Total genomic length of the candidate structure                |
| `SplitCount_sum`, `SplitCount_mean` | Junction support summary                                       |
| `CN_sum`, `CN_mean`, `CN_std`       | Copy number summary across segments                            |
| `Gene_num`                          | Number of genes overlapping the structure                      |
| `Cancergene_num`                    | Number of cancer genes overlapping the structure               |
| `SE_num`                            | Number of super-enhancer annotations overlapping the structure |

### Breakpoint-Pair File

`*_ask_breakpoint_pair.tsv` reports breakpoint pairs used during amplicon reconstruction.


| Column                      | Description                            |
| --------------------------- | -------------------------------------- |
| `Chrom1`, `Coord1`, `Clip1` | First breakpoint side and orientation  |
| `Chrom2`, `Coord2`, `Clip2` | Second breakpoint side and orientation |
| `Count`                     | Supporting read count                  |
| `Seq`                       | Junction sequence when available       |

### JCS File

`*_ask_jcs.tsv` is generated by targeted search mode.


| Column                      | Description                                                          |
| --------------------------- | -------------------------------------------------------------------- |
| `CircleID`                  | Reference circle identifier                                          |
| `total_reference_junctions` | Number of reference junctions derived from the known ecDNA structure |
| `validated_junctions`       | Number of reference junctions supported in the query BAM             |
| `total_support_reads`       | Total supporting reads across validated junctions                    |
| `JCS`                       | Junction Concordance Score                                           |
| `Detected`                  | Whether the circle passes the JCS threshold                          |

JCS is computed as:

```text
JCS = validated reference junctions / total reference junctions
```

By default, a reference junction is considered validated when it has at least five supporting reads, and a circle is marked detected when `JCS > 0.5`.

## Algorithm Overview

ASK reconstructs amplified structures from coverage and breakpoint evidence:

1. Alignment evidence extraction from an indexed BAM/CRAM.
2. Read counting in genomic bins.
3. Copy number normalization and segmentation.
4. Amplified segment detection.
5. Breakpoint detection from clipping and supplementary-alignment evidence.
6. Breakpoint-pair construction.
7. Graph-based circular and linear amplicon reconstruction.
8. Gene, cancer gene, and super-enhancer annotation.
9. ASK-style visualization.

The targeted `ask-search` workflow follows the same evidence model but constrains the initial evidence collection using a known ecDNA structure.

## Checkpointing and Modular Usage

ASK writes intermediate `.pdat` files:


| File               | Stage                                       |
| ------------------ | ------------------------------------------- |
| `*_ask_step1.pdat` | Alignment evidence and bin counts           |
| `*_ask_step2.pdat` | Copy number and amplified segment detection |
| `*_ask_step3.pdat` | Breakpoint-pair detection                   |
| `*_ask_step4.pdat` | Amplicon reconstruction                     |

These files are useful for debugging, rerunning downstream steps, and comparing parameter choices. When rerunning from scratch, use a fresh output prefix or remove incompatible intermediate files.

ASK can also be used modularly:


| Use Case                                   | Suggested Entry Point                                                 |
| ------------------------------------------ | --------------------------------------------------------------------- |
| Start from BAM/CRAM                        | `ask`                                                                 |
| Start from known ecDNA structure           | `ask-search`                                                          |
| Compare one reference ecDNA across samples | Run`ask-search` once per query BAM                                    |
| Replot existing ASK outputs                | Use the generated circular, linear, copy number, and bin-count tables |

## License

ASK is released under the MIT License. See [LICENSE](https://github.com/nanawei11/AmpliconSeeK/blob/main/LICENSE) for details. 

## Contact

For questions and feedback, please open an issue on [GitHub](https://github.com/nanawei11/AmpliconSeeK/issues) or contact Nana Wei (nnwei@shsmu.edu.cn).
