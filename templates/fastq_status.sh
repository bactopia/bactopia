#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
ERROR=0
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --extra !{extra} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{extra}
    fi
fi

if [ "!{params.skip_fastq_check}" == "false" ]; then
    # Not completely sure about the inputs, so make sure they meet minimum requirements
    echo "# fastq-scan Version" >> ${LOG_DIR}/!{task.process}.versions
    fastq-scan -v >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Check paired-end reads have same read counts
    gzip -cd !{fq[0]} | fastq-scan > r1.json
    OPTS="--sample !{sample} --min_basepairs !{params.min_basepairs} --min_reads !{params.min_reads} --min_proportion !{params.min_proportion}"
    if [ "!{single_end}" == "false" ]; then
        if ! reformat.sh in1=!{fq[0]} in2=!{fq[1]} !{qin} out=/dev/null 2> !{sample}-paired-end-error.txt; then
            ERROR=1
            echo "!{sample} FASTQs contains an error. Please check the input FASTQs.
                Further analysis is discontinued." | \
            sed 's/^\s*//' >> !{sample}-paired-end-error.txt
        else
            rm -f !{sample}-paired-end-error.txt
        fi
        gzip -cd !{fq[1]} | fastq-scan > r2.json

        if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
            ERROR=1
        fi
        rm r1.json r2.json
    else
        if ! check-fastqs.py --fq1 r1.json ${OPTS}; then
            ERROR=1
        fi
        rm r1.json
    fi
fi

if [ "${ERROR}" -eq "0" ]; then
    mkdir -p fastqs
    if [[ -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
            cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            cp  !{fq[0]} fastqs/!{sample}.fastq.gz
        fi
    fi
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
