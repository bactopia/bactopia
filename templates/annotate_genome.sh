#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p annotation
    touch annotation/!{sample}.gbk
else
    mkdir -p logs
    if [[ !{params.compress} == "true" ]]; then
        gunzip -f !{fasta}
    fi

    if [ "!{renamed}" == "true" ]; then
        echo "Original sample name (!{sample}) not used due to creating a contig ID >37 characters"
    fi

    # Prokka Version
    echo "# Prokka Version" > logs/annotate_genome.versions
    prokka --version >> logs/annotate_genome.versions 2>&1
    prokka --outdir annotation \
        --force \
        --prefix '!{sample}' \
        --genus '!{genus}' \
        --species '!{species}' \
        --evalue '!{params.prokka_evalue}' \
        --coverage !{params.prokka_coverage} \
        --cpus !{task.cpus} \
        --centre '!{params.centre}' \
        --mincontiglen !{params.min_contig_len} \
        !{locustag} \
        !{prodigal} \
        !{addgenes} \
        !{compliant} \
        !{proteins} \
        !{rawproduct} \
        !{cdsrnaolap} \
        !{addmrna} \
        !{norrna} \
        !{notrna} \
        !{rnammer} \
        !{rfam} \
        !{gunzip_fasta} > logs/prokka.out 2> logs/prokka.err

    if [[ !{params.compress} == "true" ]]; then
        find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err logs/annotate_genome.err
        cp .command.out logs/annotate_genome.out
        cp .command.run logs/annotate_genome.run
        cp .command.sh logs/annotate_genome.sh
        cp .command.trace logs/annotate_genome.trace
    else
        rm -rf logs/
    fi
fi
