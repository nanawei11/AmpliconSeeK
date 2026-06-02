#!/bin/bash

# set path of the pipeline
bdir=`dirname $0`

# input
downloaded=$1    
mapped=${downloaded/_downloaded.tsv/_mapped.tsv} 

# loop mapping
while read line; do
  
    item=($line)   # ERP116362_ERR3452702_BRCA_Hs578T_Input SRRid
    bname=${item[0]}    # ERP116362_ERR3452702_BRCA_Hs578T_Input
    log_file=${bname}/log/map_bwa.log   # ERP116362_ERR3452702_BRCA_Hs578T_Input/log/map_bwa.log

    # create log file
    if [ ! -d ${bname}/log ]; then
        mkdir -p ${bname}/log
    fi

    # add time stamp
    echo -e "\nProcess date: `date`\n" >> $log_file

    # run main
    $bdir/map_bwa.sh $bname $mapped &>> $log_file

done < $downloaded


