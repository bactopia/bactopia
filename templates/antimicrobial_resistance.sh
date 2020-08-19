#!/bin/bash
set -e
set -u

LOG_DIR="logs/!{task.process}"
mkdir -p ${LOG_DIR}
touch ${LOG_DIR}/!{task.process}.versions

if [[ !{params.compress} == "true" ]]; then
    zcat !{genes} > !{sample}.ffn
    zcat !{proteins} > !{sample}.faa
fi

tar -xzvf !{amrdb}
mkdir antimicrobial_resistance

# amrfinder Version
echo "# amrfinder Version" >> ${LOG_DIR}/!{task.process}.versions
amrfinder --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
amrfinder -n !{sample}.ffn \
          -d amrdb/ \
          -o antimicrobial_resistance/!{sample}-gene-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{task.cpus} !{organism_gene} !{plus} !{report_common} > ${LOG_DIR}/amrfinder-gene.out 2> ${LOG_DIR}/amrfinder-gene.err

amrfinder -p !{sample}.faa \
          -d amrdb/ \
          -o antimicrobial_resistance/!{sample}-protein-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{task.cpus} !{organism_protein} !{plus} !{report_common} > ${LOG_DIR}/amrfinder-protein.out 2> ${LOG_DIR}/amrfinder-protein.err

if [[ !{params.compress} == "true" ]]; then
    rm !{sample}.faa !{sample}.ffn
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.run ${LOG_DIR}/!{task.process}.run
    cp .command.sh ${LOG_DIR}/!{task.process}.sh
    cp .command.trace ${LOG_DIR}/!{task.process}.trace
else
    rm -rf ${LOG_DIR}/
fi
