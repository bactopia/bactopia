#!/bin/bash
set -e
set -u

shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize !{params.genome_size} --outdir . \
    --minlen 500 --cpus !{cpus} --assembler !{params.assembler} --noreadcorr --force
mv contigs.fa !{sample}.fna
assembly-scan !{sample}.fna > !{sample}.fna.json
gzip !{sample}.fna
