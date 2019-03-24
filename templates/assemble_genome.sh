#!/bin/bash
set -e
set -u

GENOME_SIZE=`head -n 1 !{genome_size_file}`

if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize ${GENOME_SIZE} --outdir . \
        --minlen 500 --cpus !{cpus} --assembler !{params.assembler} --noreadcorr --force
else
    # Single-End Reads
    shovill-se --se !{fq[0]} --depth 0 --gsize ${GENOME_SIZE} --outdir . \
        --minlen 500 --cpus !{cpus} --assembler !{params.assembler} --force
fi

mv contigs.fa !{sample}.fna
assembly-scan !{sample}.fna > !{sample}.fna.json
gzip !{sample}.fna
