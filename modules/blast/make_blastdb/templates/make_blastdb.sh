#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir blastdb
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# makeblastdb Version" >> ${LOG_DIR}/!{task.process}.versions
makeblastdb -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

# Verify AWS files were staged
if [[ ! -L "!{fasta}" ]]; then
    check-staging.py --assembly !{fasta}
fi

if [[ !{params.compress} == "true" ]]; then
    gzip -cd !{fasta} | \
    makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
else
    cat !{fasta} | \
    makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
