#!/bin/bash

# set path of the pipeline
bdir="$(dirname "$0")"

# Check input arguments
if [ $# -lt 1 ]; then
    echo "Usage: bash $0 <srr_table.tsv>"
    echo "Example: bash $0 metadata.tsv"
    exit 1
fi

# input
srrtable="$1"
downloaded="${srrtable/.tsv/_downloaded.tsv}"

if [ ! -e "$srrtable" ]; then
    echo "[ERROR] Cannot find SRR table: $srrtable"
    exit 1
fi

# loop of SRR to fastq
while read -r sid srrs extra; do

    # Skip empty lines
    if [ -z "${sid:-}" ]; then
        continue
    fi

    # Skip comment lines
    if [[ "$sid" =~ ^# ]]; then
        continue
    fi

    # Skip possible header line
    if [[ "$sid" == "sid" || "$sid" == "sample" || "$sid" == "sample_id" ]]; then
        continue
    fi

    # Check SRR column
    if [ -z "${srrs:-}" ]; then
        echo "[WARN] Missing SRR ID for sample: $sid"
        continue
    fi

    # Reset failed flag for each sample
    download_failed=0

    # directories
    fdir="${sid}/fastq"
    ldir="${sid}/log"
    log_file="${ldir}/fastq_download.log"

    # create run folder and log folder
    mkdir -p "$fdir"
    mkdir -p "$ldir"

    # create a file recording conflict run
    if [ ! -e "${fdir}/conflict.txt" ]; then
        touch "${fdir}/conflict.txt"
    fi

    {
        echo
        echo "============================================================"
        echo "[INFO] Process date: $(date)"
        echo "[INFO] Sample ID: $sid"
        echo "[INFO] SRR list: $srrs"
        echo "[INFO] FASTQ dir: $fdir"
        echo "============================================================"
        echo
    } >> "$log_file"

    # download SRRs
    IFS='|' read -ra srr_array <<< "$srrs"

    for srr in "${srr_array[@]}"; do

        # Skip empty SRR entries
        if [ -z "$srr" ]; then
            continue
        fi

        echo "[INFO] Downloading SRR: $srr"

        {
            echo
            echo "------------------------------"
            echo "[INFO] Start SRR: $srr"
            echo "[INFO] Date: $(date)"
            echo "------------------------------"
        } >> "$log_file"

        bash "${bdir}/srr2fq.sh" "$srr" "$sid" >> "$log_file" 2>&1

        if [ $? -ne 0 ]; then
            echo "[ERROR] Failed to download SRR: $srr"
            echo "[ERROR] Failed to download SRR: $srr" >> "$log_file"
            download_failed=1
        else
            echo "[INFO] Finished SRR: $srr"
            echo "[INFO] Finished SRR: $srr" >> "$log_file"
        fi

    done

    # run check / merge
    if [ "$download_failed" -eq 0 ]; then
        echo "[INFO] Running fastqMerge.sh for sample: $sid"
        echo "[INFO] Running fastqMerge.sh for sample: $sid" >> "$log_file"

        bash "${bdir}/fastqMerge.sh" "$srrs" "$sid" "$downloaded" >> "$log_file" 2>&1

        if [ $? -ne 0 ]; then
            echo "[ERROR] fastqMerge.sh failed for sample: $sid"
            echo "[ERROR] fastqMerge.sh failed for sample: $sid" >> "$log_file"
        else
            echo "[INFO] fastqMerge.sh finished for sample: $sid"
            echo "[INFO] fastqMerge.sh finished for sample: $sid" >> "$log_file"
        fi
    else
        echo "[WARN] Skip fastqMerge.sh because at least one SRR failed: $sid"
        echo "[WARN] Skip fastqMerge.sh because at least one SRR failed: $sid" >> "$log_file"
    fi

    # remove the conflict record
    rm -f "${fdir}/conflict.txt"

done < "$srrtable"