#!/bin/bash

# set path of the pipeline
bdir=`dirname $0`

# input
ask=/cluster2/home/WeiNa/project/ecDNA/ask_run/ask3/ask_cmd.py
bname=$1
signal_file=$2

## pars
bam=$bname/bwa/${bname}.bam
bw=$bname/bwa/${bname}.bw
askdir=ask
outprefix=$bname/$askdir/$bname

## make dir for ask
if [ ! -d $bname/$askdir ]
then
    mkdir $bname/$askdir
fi

# create a file recording conflict run

if [ -e $bname/$askdir/conflict.txt ]; then
    exit 1
else
    touch $bname/$askdir/conflict.txt
fi


t=0
if [ -e ${outprefix}_ask_step1.pdat ]
then
    t=1
fi

if [  -e ${outprefix}_ask_step2.pdat ]
then
    t=2
fi

if [  -e ${outprefix}_ask_step3.pdat ]
then
    t=3
fi

# # run ask 
# if [ ! -e ${outprefix}_ask_stats.tsv ]
# then
#     $ask -i $bam -o $outprefix --knn 3 --subseg --method both --rm_amp_df --juncread 5  --SA_with_nm
# fi

if [ ! -e ${outprefix}_ask_stats.tsv ]
then
    # if [ -e ${outprefix}_ask_step3.pdat ]
    lowercase=$(echo "$bam" | tr '[:upper:]' '[:lower:]')
    if [[ $lowercase = *"input"*  ]] 
    then
        echo $bam stand
        ~/anaconda3/envs/ask/bin/python $ask -i $bam -o $outprefix --knn 3 --subseg --method both --juncread 5  --SA_with_nm --run_from_pdat ${t}
    else
        echo $bam bais
        ~/anaconda3/envs/ask/bin/python $ask -i $bam -o $outprefix --knn 3 --subseg --method both --juncread 5 -d bias --SA_with_nm --run_from_pdat ${t}
   
    fi 
fi

## raise signal
if [ -e ${outprefix}_ask_stats.tsv ]; then
    if [ -e $signal_file ]; then
        if [ -z `awk '$0 == "'${bname}'"' $signal_file` ]; then
            echo "${bname}" >> $signal_file
        fi
    else
        echo "${bname}" > $signal_file
    fi
fi

# remove a file recording conflict run
if [ -e $bname/$askdir/conflict.txt ]; then
    # remove the conflict record
    rm $bname/$askdir/conflict.txt
fi
