#!/bin/bash
set -e
set -u

snippy !{fastq} --ref !{reference} --prefix !{sample} \
    --outdir !{reference_name} --cpus !{task.cpus}

vcf-annotator !{reference_name}/!{sample}.vcf !{reference} > !{reference_name}/!{sample}-final.vcf

rm -rf !{reference_name}/reference !{reference_name}/ref.fa* !{reference_name}/!{sample}.vcf.gz*

find !{reference_name}/ -type f -not -name "*.bam*" -and -not -name "*.log*" -and -not -name "*.txt*" | \
    xargs -I {} pigz -n --best -p !{task.cpus} {}
