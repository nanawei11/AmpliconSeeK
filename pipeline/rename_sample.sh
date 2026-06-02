#!/bin/bash
#####
# rename samples based on geo id
# used to reformat sample names
#####

# set path of the pipeline
bdir=`dirname $0`

# input
srrtable=$1

# loop of SRR to fastq
while read line; do
    
    item=($line)
    sid=${item[0]}
    bb=(${sid//_/ })
    bname=${bb[0]}_${bb[1]}
    sid_old=`ls --color=no |grep $bname`

    # change the root folder
    mv $sid_old $sid
    cd $sid

    # change name of files involved
    paste <(find . -name "$bname*") \
        <(find . -name "$bname*" \
        |sed "s/$sid_old/$sid/g") \
        |sed 's/\t/ /g' |sed 's/^/mv /g' |bash

    # cd back
    cd -
done < $srrtable


