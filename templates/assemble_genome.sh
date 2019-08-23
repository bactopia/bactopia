#!/bin/bash
set -e
set -u

GENOME_SIZE=`head -n 1 !{genome_size}`

if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize ${GENOME_SIZE} \
        --outdir . \
        --force \
        --minlen !{params.min_contig_len} \
        --mincov !{params.min_contig_cov} \
        --namefmt "!{params.contig_namefmt}" \
        --keepfiles \
        --cpus !{task.cpus} \
        --ram !{shovill_ram} \
        --assembler !{params.assembler} \
        --noreadcorr !{opts} !{kmers} !{nostitch} !{nocorr}
else
    # Single-End Reads
    shovill-se --se !{fq[0]} --depth 0 --gsize ${GENOME_SIZE} \
        --outdir . \
        --force \
        --minlen !{params.min_contig_len} \
        --mincov !{params.min_contig_cov} \
        --namefmt "!{params.contig_namefmt}" \
        --keepfiles \
        --cpus !{task.cpus} \
        --ram !{shovill_ram} \
        --assembler !{params.assembler} !{opts} !{kmers} !{nocorr}
fi

TOTAL_CONTIGS=`grep -c "^>" contigs.fa || true`
if [ "${TOTAL_CONTIGS}" -gt "0" ]; then
    mv contigs.fa !{sample}.fna
    assembly-scan !{sample}.fna > !{sample}.fna.json
    gzip !{sample}.fna
else
    echo "Assembly was successful, but 0 contigs were formed." > shovill-zero-contigs.txt
fi

if [ "!{params.keep_all_files}" == "false" ]; then
    # Remove intermediate files
    rm -fv shovill.bam* flash.extendedFrags* flash.notCombined* skesa.fasta.*
fi
