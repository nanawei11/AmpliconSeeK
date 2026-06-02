#!/bin/bash

# set path of the pipeline
bdir=`dirname $0`

# input
ask=ask_cmd.py
bname=$1
signal_file=$2
opt=$3

## pars
bam=$bname/bwa/$bname.bam
bw=$bname/bwa/$bname.bw
outprefix=$bname/ask/$bname

## make dir for ask
if [ ! -d $bname/ask ]
then
    mkdir $bname/ask
fi

# create a file recording conflict run
if [ -e $bname/ask/conflict.txt ]; then
    exit 1
else
    touch $bname/ask/conflict.txt
fi

## run ask
# if [ ! -e ${outprefix}_ask_stats.csv ]
# then
# python $ask -i $bam -b $bw -o $outprefix $opt
# fi
python $ask -i $bam -b $bw -o $outprefix $opt

## raise signal
if [ -e $signal_file ]; then
    if [ -z `awk '$0 == "'${bname}'"' $signal_file` ]; then
        echo "${bname}" >> $signal_file
    fi
else
	echo "${bname}" > $signal_file
fi

# remove the conflict record
rm $bname/ask/conflict.txt

