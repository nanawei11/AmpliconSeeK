#!/bin/bash

# set path of the pipeline
bdir=`dirname $0`

# input
mapped=$1
signal_file=${mapped/_mapped.tsv/_asked.tsv}      # sample/retrun_Input_asked.tsv

# loop of SRR to fastq
while read line; do
  
    item=($line)
    bname=${item[0]}
    askdir=$bname/ask
    log_file=${bname}/log/ask

    if [ -e $bname/ask/conflict.txt ]; then
        continue
    fi

    if [ ! -e $signal_file ]; then
        touch $signal_file
    fi

    # raise signal
    if [ -e $signal_file ]; then
        if [ -z `awk '$0 == "'${bname}'"' $signal_file` ]; then
            if [ -d $askdir ]; then
                # rm -rf $askdir
                rm -rf $log_file
                echo "Sample: <$bname> - re-asking..."
            fi
        fi
    fi

    # create log file
    if [ ! -d ${bname}/log ]; then
        mkdir -p ${bname}/log
    fi

    # add time stamp
    echo -e "\nProcess date: `date`\n" >> $log_file

    ## run ask
    echo $bname $signal_file
    sh $bdir/ask.sh $bname $signal_file &>> $log_file

done < $mapped


