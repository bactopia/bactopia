#!/bin/bash
set -e
set -u

snippy !{fastq} \
    --ref !{reference} \
    --cpus !{task.cpus} \
    --ram !{snippy_ram} \
    --outdir !{reference_name} \
    --prefix !{sample} \
    --mapqual !{params.mapqual} \
    --basequal !{params.basequal} \
    --mincov !{params.mincov} \
    --minfrac !{params.minfrac} \
    --minqual !{params.minqual} \
    --maxsoft !{params.maxsoft} !{bwaopt} !{fbopt}

# Add GenBank annotations
vcf-annotator !{reference_name}/!{sample}.vcf !{reference} > !{reference_name}/!{sample}.annotated.vcf

# Get per-base coverage
grep "^##contig" !{reference_name}/!{sample}.vcf > !{sample}.coverage.txt
genomeCoverageBed -ibam !{reference_name}/!{sample}.bam -d | cut -f3 >> !{reference_name}/!{sample}.coverage.txt

# Mask low coverage regions
mask-consensus.py !{sample} !{reference_name} \
                  !{reference_name}/!{sample}.consensus.subs.fa \
                  !{reference_name}/!{sample}.subs.vcf \
                  !{reference_name}/!{sample}.coverage.txt > !{reference_name}/!{sample}.masked-consensus.subs.fa

# Clean Up
rm -rf !{reference_name}/reference !{reference_name}/ref.fa* !{reference_name}/!{sample}.vcf.gz*

if [[ !{params.compress} == "true" ]]; then
    find !{reference_name}/ -type f -not -name "*.bam*" -and -not -name "*.log*" -and -not -name "*.txt*" | \
        xargs -I {} pigz -n --best -p !{task.cpus} {}
    pigz -n --best -p !{task.cpus} !{reference_name}/!{sample}.coverage.txt
fi
