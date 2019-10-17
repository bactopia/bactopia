#!/bin/bash
set -e
set -u

if [[ !{params.compress} == "true" ]]; then
    zcat !{genes} > !{sample}.ffn
    zcat !{proteins} > !{sample}.faa
fi

mkdir antimicrobial_resistance
amrfinder -n !{sample}.ffn \
          -o antimicrobial_resistance/!{sample}-gene-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{task.cpus} !{organism_gene} !{plus} !{report_common}

amrfinder -p !{sample}.faa \
          -o antimicrobial_resistance/!{sample}-protein-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{task.cpus} !{organism_protein} !{plus} !{report_common}

if [[ !{params.compress} == "true" ]]; then
    rm !{sample}.faa !{sample}.ffn
fi

