#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}/!{dataset_name}"
if [ "!{params.dry_run}" == "true" ]; then
    touch minmer_query.!{dataset_name}.dry_run.txt
else
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

    if [ "!{dataset_name}" == "refseq-k21-s1000.msh" ]; then
        echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
        mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

        printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-refseq-k21.txt
        zcat !{fastq} | \
        mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
        sort -gr >> !{sample}-refseq-k21.txt 2> ${LOG_DIR}/mash.err
    elif [ "!{dataset_name}" == "plsdb.msh" ]; then
        echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
        mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

        printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-plsdb-k21.txt
        zcat !{fastq} | \
        mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
        sort -gr >> !{sample}-plsdb-k21.txt 2> ${LOG_DIR}/mash.err
    elif [ "!{dataset_name}" == "genbank-k21.json.gz" ]; then
        echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
        sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
        sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k21.txt 2> ${LOG_DIR}/sourmash.err
    elif [ "!{dataset_name}" == "genbank-k31.json.gz" ]; then
        echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
        sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
        sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k31.txt 2> ${LOG_DIR}/sourmash.err
    else
        echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
        sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
        sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k51.txt 2> ${LOG_DIR}/sourmash.err
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
fi
