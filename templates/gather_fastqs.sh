#!/bin/bash
set -e
set -u

mkdir fastqs
if [ "!{single_end}" == "is_accession" ]; then
    # Download accession from ENA
    ena-dl !{sample} fastqs/ --group_by_experiment --is_experiment
else
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
    fi
fi
