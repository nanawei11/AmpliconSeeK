#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: bash $0 <downloaded.tsv>"
    echo "Example: bash $0 metadata_downloaded.tsv"
    exit 1
fi

# input
bname=$1   # ERP116362_ERR3452702_BRCA_Hs578T_Input
signal_file=$2  # sample/rerun_mapped.tsv

## pars
nCPU=10 

bwa_index=/cluster2/home/WeiNa/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa 
# picard_cmd=/michorlab/hjwu/soft/picard-2.22.4/picard.jar
fastq=`ls -d $bname/fastq/*gz |xargs`  # get the fastq file
bam=$bname/bwa/${bname}.bam        # ERP116362_ERR3452702_BRCA_Hs578T_Input/bwa/ERP116362_ERR3452702_BRCA_Hs578T_Input.bam
bam_unsorted=$bname/bwa/$bname.presort.bam           # ERP116362_ERR3452702_BRCA_Hs578T_Input/bwa/ERP116362_ERR3452702_BRCA_Hs578T_Input.presort.bam
tmp=$bname/bwa/$bname             # ERP116362_ERR3452702_BRCA_Hs578T_Input/bwa/ERP116362_ERR3452702_BRCA_Hs578T_Input
bigwig=$bname/bwa/${bname}.bw       # ERP116362_ERR3452702_BRCA_Hs578T_Input/bwa/ERP116362_ERR3452702_BRCA_Hs578T_Input.bw
bambai=$bname/bwa/${bname}.bam.bai
samtools=~/anaconda3/envs/ask/bin/samtools
## make dir for bwa:  ERP116362_ERR3452702_BRCA_Hs578T_Input/bwa
if [ ! -d $bname/bwa ]
then
    mkdir $bname/bwa
fi

# create a file recording conflict run
if [ -e $bname/bwa/conflict.txt ]; then
    exit 1
else
    touch $bname/bwa/conflict.txt
fi

## run mapping
if [ ! -e $bam ] && [ ! -e $bam_unsorted ]
then
    bwa mem -t $nCPU $bwa_index $fastq | \
        $samtools view -Shb - > $bam_unsorted
fi

## run sorting
if [ ! -e $bam ]
then
    $samtools fixmate --threads $nCPU -m $bam_unsorted - \
        |$samtools sort --threads $nCPU -T $tmp - \
        |$samtools markdup --threads $nCPU -T $tmp -S -s - $bam
    
    $samtools index $bam
    bamCoverage -p $nCPU --ignoreDuplicates -b $bam -o $bigwig
fi

## raise signal
if [ -e $signal_file ]; then
    if [ -z `awk '$0 == "'${bname}'"' $signal_file` ]; then  # STRING 的长度为零则为真
        echo "${bname}" >> $signal_file
    fi
else
    echo "${bname}" > $signal_file
fi

## remove tmp file
rm $bam_unsorted

# remove the conflict record
rm $bname/bwa/conflict.txt



