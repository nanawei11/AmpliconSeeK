#####
# batch merge multiple fastq.gz files into one file.
#####

# set path of the pipeline
bdir=`dirname $0`  

# input
srrs=$1  # SRRXXX
sid=$2   # ERP116362_ERR3452702_BRCA_Hs578T_Input
signal_file=$3  # sample/rerun_downloaded.tsv
fdir=${sid}/fastq   # ERP116362_ERR3452702_BRCA_Hs578T_Input/fastq
outprefix=${sid}

# set pars
x=0
y=0
z=0
R1=""
R2=""

# run
if [ ! -e $fdir/${outprefix}_R1.fastq.gz ]; then

    # check for all SRR files of a GSM are all downloaded
    for srr in $(echo $srrs | tr "|" "\n"); do
        fdir=${sid}/fastq
        ## read R1
        if [ -e $fdir/${srr}_R1.fastq.gz ]; then
            JJ1="Correct"
            R1+=" $fdir/${srr}_R1.fastq.gz"
            let "y+=1"
        else
            JJ1="Error"
            x=1
        fi
        ## read R2
        if [ -e $fdir/${srr}_R2.fastq.gz ]; then
            JJ2="Correct"
            R2+=" $fdir/${srr}_R2.fastq.gz"
            let "z+=1"
        else
            JJ2="Error"
        fi
        echo -e $srr"\t"$JJ1"\t"$JJ2
    done

    # echo "$x $y $z $fdir/${sid}_R1.fastq.gz"
    # if all files are downloaded successfully, start to merge fastq files
    if [ $x == 0 ]; then

        if [ $y == 1 ]; then
            mv $R1 $fdir/${sid}_R1.fastq.gz
        elif [ $y -gt 1 ]; then
            cat $R1 > $fdir/${sid}_R1.fastq.gz
            rm $R1
        fi

        if [ $z == 1 ]; then
            mv $R2 $fdir/${sid}_R2.fastq.gz
        elif [ $z -gt 1 ]; then
            cat $R2 > $fdir/${sid}_R2.fastq.gz
            rm $R2
        fi

        ## raise signal only if the fastq files are successfully merged
        if [ -e $signal_file ]; then
            if [ -z `awk '$0 == "'${sid}'"' $signal_file` ]; then
                echo "${sid}" >> $signal_file
            fi
        else
            echo "${sid}" > $signal_file
        fi

        echo ">>> SRR files are downloaded and merged successfully!"
    else
        echo ">>> Some of the SRR files are NOT downloaded successfully!"
    fi
else
    ## raise signal
    echo "/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq/pipeline/fastqMerge.sh"
    echo $signal_file
    echo ${sid}
    if [ -e $signal_file ]; then
        if [ -z `awk '$0 == "'${sid}'"' $signal_file` ]; then
            echo "${sid}" >> $signal_file
        fi
    else
        echo ${sid}
        echo "${sid}" > $signal_file
    fi
fi


