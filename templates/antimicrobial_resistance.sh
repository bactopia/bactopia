#!/bin/bash
set -e
set -u

mkdir antimicrobial_resistance
zcat !{genes} > !{sample}.ffn
amrfinder -n !{sample}.ffn \
          -o antimicrobial_resistance/!{sample}-gene-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{cpus} !{organism_gene} !{plus} !{report_common}

zcat !{proteins} > !{sample}.faa
amrfinder -p !{sample}.faa \
          -o antimicrobial_resistance/!{sample}-protein-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{cpus} !{organism_protein} !{plus} !{report_common}

rm !{sample}.faa !{sample}.ffn
