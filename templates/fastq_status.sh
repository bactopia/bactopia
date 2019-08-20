#!/bin/bash
set -e
set -u

SEQUENCED_BP=`zcat *.gz | fastq-scan | grep "total_bp" | sed -r 's/.*:([0-9]+),/\1/'`
TOTAL_READS=`zcat *.gz | fastq-scan | grep "read_total" | sed -r 's/.*:([0-9]+),/\1/'`
if [ ${SEQUENCED_BP} -gt "!{params.min_basepairs}" ] && [ ${TOTAL_READS} -gt "!{params.min_reads}" ]; then
    mkdir -p fastqs
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
    fi
else
    if [ ${SEQUENCED_BP} -lt "!{params.min_basepairs}" ]; then
        mkdir -p quality-control
        echo "!{sample} FASTQ(s) contain ${SEQUENCED_BP} total basepairs. This does not
              exceed the required minimum !{params.min_basepairs} bp. Further analysis is
              discontinued." | \
        sed 's/^\s*//' > low-sequence-depth-error.txt
    fi

    if [ ${TOTAL_READS} -lt "!{params.min_reads}" ]; then
        mkdir -p quality-control
        echo "!{sample} FASTQ(s) contain ${TOTOAL_READS} total basepairs. This does not
              exceed the required minimum !{params.min_reads} read count. Further analysis is
              discontinued." | \
        sed 's/^\s*//' > low-read-count-error.txt
    fi
fi
