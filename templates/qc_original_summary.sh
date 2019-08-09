#!/bin/bash
set -e
set -u

GENOME_SIZE=`head -n 1 !{genome_size}`
if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    zcat !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R1-original.json
    zcat !{fq[1]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R2-original.json
    ln -s !{fq[0]} !{sample}_R1-original.fastq.gz
    ln -s !{fq[1]} !{sample}_R2-original.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}_R1-original.fastq.gz !{sample}_R2-original.fastq.gz
else
    # Single-End Reads
    zcat !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}-original.json
    ln -s !{fq[0]} !{sample}-original.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}-original.fastq.gz
fi

mkdir -p quality-control/summary-original
mv *.json  quality-control/summary-original
mv *fastqc.html quality-control/summary-original
mv *fastqc.zip quality-control/summary-original
