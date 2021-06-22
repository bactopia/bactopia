#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}/

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
if [[ !{params.compress} == "true" ]]; then
    gunzip -f !{fasta}
fi

if [ "!{renamed}" == "true" ]; then
    echo "Original sample name (!{sample}) not used due to creating a contig ID >37 characters"
fi

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --assembly !{gunzip_fasta} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --assembly !{gunzip_fasta}
    fi
fi

# Prokka Version
echo "# Prokka Version" >> ${LOG_DIR}/!{task.process}.versions
prokka --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
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
    !{gunzip_fasta} > ${LOG_DIR}/prokka.out 2> ${LOG_DIR}/prokka.err

if [[ !{params.compress} == "true" ]]; then
    find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
        xargs -I {} pigz -n --best -p !{task.cpus} {}
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
