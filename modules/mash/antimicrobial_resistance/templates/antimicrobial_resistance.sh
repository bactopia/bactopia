#!/bin/bash
set -e
set -u

LOG_DIR="logs/!{task.process}"
mkdir -p ${LOG_DIR}

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Verify AWS files were staged
if [[ ! -L "!{genes} " ]]; then
    check-staging.py --fq1 !{genes} --fq2 !{proteins} --extra !{amrdb}
fi

if [[ !{params.compress} == "true" ]]; then
    gzip -cd !{genes} > !{sample}.ffn
    gzip -cd !{proteins} > !{sample}.faa
fi

tar -xzvf !{amrdb}
mkdir !{amrdir}

# amrfinder Version
echo "# amrfinder Version" >> ${LOG_DIR}/!{task.process}.versions
amrfinder --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
amrfinder -n !{sample}.ffn \
          -d amrfinderdb/ \
          -o !{amrdir}/!{sample}-gene-report.txt \
          --ident_min !{params.amr_ident_min} \
          --coverage_min !{params.amr_coverage_min} \
          --translation_table !{params.amr_translation_table} \
          --threads !{task.cpus} !{organism_gene} !{plus} !{report_common} > ${LOG_DIR}/amrfinder-gene.out 2> ${LOG_DIR}/amrfinder-gene.err

amrfinder -p !{sample}.faa \
          -d amrfinderdb/ \
          -o !{amrdir}/!{sample}-protein-report.txt \
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
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
