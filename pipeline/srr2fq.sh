#!/bin/bash

# Check input arguments
if [ $# -lt 2 ]; then
    echo "Usage: sh $0 <SRR_ID> <SAMPLE_ID>"
    echo "Example: sh $0 SRR5661617 example_sample"
    exit 1
fi

# input
srr=$1   # SRR ID
sid=$2   # Sample name
fdir=${sid}/fastq
outprefix=${sid}

if [ ! -d $fdir ]; then
    mkdir -p $fdir
fi

if [ ! -e $fdir/${outprefix}_R1.fastq.gz ]; then
    
    if [ ! -e $fdir/${srr}_R1.fastq.gz ]; then
        
        if [ ! -e $fdir/${srr}_1.fastq.gz ]; then
                

            echo "[INFO] Trying parallel-fastq-dump for $srr"

            parallel-fastq-dump -t 12 --split-files --gzip -O $fdir -s $srr 

            
            if [ $? -ne 0 ]; then
                echo "[WARN] parallel-fastq-dump failed for $srr"
                echo "[INFO] Falling back to prefetch + fastq-dump"


                prefetch --max-size 2000000000000 $srr

                if [ $? -ne 0 ]; then
                    echo "[ERROR] prefetch failed for $srr"
                    exit 1
                fi


                vdb-validate $srr

                if [ $? -ne 0 ]; then
                    echo "[ERROR] vdb-validate failed for $srr"
                    exit 1
                fi

                fastq-dump --gzip --split-files -O $fdir $srr

                if [ $? -ne 0 ]; then
                    echo "[ERROR] fastq-dump failed for $srr"
                    exit 1
                fi


                # clean
                rm -rf $srr

            fi


            if [ -e $fdir/${srr}_1.fastq.gz ]; then
                mv $fdir/${srr}_1.fastq.gz $fdir/${srr}_R1.fastq.gz
            fi

            if [ -e $fdir/${srr}_2.fastq.gz ]; then
                mv $fdir/${srr}_2.fastq.gz $fdir/${srr}_R2.fastq.gz
            fi

        fi
            
    fi

fi