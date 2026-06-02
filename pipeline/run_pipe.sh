#!/bin/bash
#####
# run pipe
#####

# set path of the pipeline
bdir=`dirname $0`  # pipeline/

## input
srrtable=$1
downloaded=${srrtable/.tsv/_downloaded.tsv}  
mapped=${srrtable/.tsv/_mapped.tsv}
# mapped=$1
asked=${srrtable/.tsv/_asked.tsv}

## download data
$bdir/srr2fq_batch.sh $srrtable &> $bdir/log/log_srr2fq 

## mapping
sh $bdir/map_bwa_batch.sh $downloaded &> $bdir/log/log_map_bwa &  # pipeline/log/

## ask
sh $bdir/ask_batch.sh $mapped #&> $bdir/log/log_ask    record the ask_batch.sh 
