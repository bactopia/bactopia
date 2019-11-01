#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p fastqs
    touch fastqs/!{sample}.fastq.gz
else
    ERROR=0
    SEQUENCED_BP=`zcat *.gz | fastq-scan | grep "total_bp" | sed -r 's/.*:([0-9]+),/\1/'`
    TOTAL_READS=`zcat *.gz | fastq-scan | grep "read_total" | sed -r 's/.*:([0-9]+),/\1/'`

    # Check paired-end reads have same read counts
    if [ "!{single_end}" == "false" ]; then
        R1_COUNT=`zcat !{fq[0]} | wc -l`
        R2_COUNT=`zcat !{fq[1]} | wc -l`
        if [ "${R1_COUNT}" != "${R2_COUNT}" ]; then
            ERROR=1
            echo "!{sample} FASTQs contain different read counts. Please check the input FASTQs.
                  Further analysis is discontinued." | \
            sed 's/^\s*//' > different-read-counts-error.txt
        fi
    fi

    if [ ${SEQUENCED_BP} -lt "!{params.min_basepairs}" ]; then
        ERROR=1
        echo "!{sample} FASTQ(s) contain ${SEQUENCED_BP} total basepairs. This does not
              exceed the required minimum !{params.min_basepairs} bp. Further analysis is
              discontinued." | \
        sed 's/^\s*//' > low-sequence-depth-error.txt
    fi

    if [ ${TOTAL_READS} -lt "!{params.min_reads}" ]; then
        echo "!{sample} FASTQ(s) contain ${TOTOAL_READS} total basepairs. This does not
              exceed the required minimum !{params.min_reads} read count. Further analysis is
              discontinued." | \
        sed 's/^\s*//' > low-read-count-error.txt
    fi

    if [ "${ERROR}" == "0" ]; then
        mkdir -p fastqs
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    fi
fi
