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

vcf-annotator !{reference_name}/!{sample}.vcf !{reference} > !{reference_name}/!{sample}.annotated.vcf

genomeCoverageBed -ibam !{reference_name}/!{sample}.bam -d > !{reference_name}/!{sample}.coverage.txt

# Clean Up
rm -rf !{reference_name}/reference !{reference_name}/ref.fa* !{reference_name}/!{sample}.vcf.gz*

if [[ !{params.compress} == "true" ]]; then
    find !{reference_name}/ -type f -not -name "*.bam*" -and -not -name "*.log*" -and -not -name "*.txt*" | \
        xargs -I {} pigz -n --best -p !{task.cpus} {}
    pigz -n --best -p !{task.cpus} !{reference_name}/!{sample}.coverage.txt
fi
