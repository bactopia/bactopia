#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p annotation
    touch annotation/!{sample}.gbk
else
    if [[ !{params.compress} == "true" ]]; then
        gunzip -f !{fasta}
    fi

    prokka --outdir annotation \
        --force \
        --prefix !{sample} \
        !{addgenes} \
        !{addmrna} \
        --locustag !{sample} \
        --centre !{params.centre} \
        --genus !{genus} \
        --species !{species} \
        !{proteins} \
        !{rawproduct} \
        !{cdsrnaolap} \
        --evalue '!{params.prokka_evalue}' \
        --coverage !{params.prokka_coverage} \
        --cpus !{task.cpus} \
        --mincontiglen !{params.min_contig_len} \
        !{norrna} \
        !{notrna} \
        !{rnammer} \
        !{rfam} \
        !{gunzip_fasta}

    if [[ !{params.compress} == "true" ]]; then
        find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
    fi
fi
