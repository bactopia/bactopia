#!/bin/bash
set -e
set -u

gunzip -f !{fasta}
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
    !{gunzip_fasta}

find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
    xargs -I {} pigz -n -p !{task.cpus} {}
