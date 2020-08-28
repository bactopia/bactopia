#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p fastqs
    touch fastqs/!{sample}.fastq.gz
else
    ERROR=0
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    if [ "!{params.skip_fastq_check}" == "false" ]; then
        # Not completely sure about the inputs, so make sure they meet minimum requirements
        echo "# fastq-scan Version" >> ${LOG_DIR}/!{task.process}.versions
        fastq-scan -v >> ${LOG_DIR}/!{task.process}.versions 2>&1
        gzip -cd *.fastq.gz | fastq-scan > info.txt
        SEQUENCED_BP=`grep "total_bp" info.txt | sed -r 's/.*: ([0-9]+),/\1/'`
        TOTAL_READS=`grep "read_total" info.txt | sed -r 's/.*: ([0-9]+),/\1/'`
        rm info.txt

        # Check paired-end reads have same read counts
        if [ "!{single_end}" == "false" ]; then
            if ! reformat.sh in1=!{fq[0]} in2=!{fq[1]} !{qin} out=/dev/null 2> !{sample}-paired-end-error.txt; then
                ERROR=1
                echo "!{sample} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \
                sed 's/^\s*//' >> !{sample}-paired-end-error.txt
            else
                rm -f !{sample}-paired-end-error.txt
            fi
        fi

        if [ ${SEQUENCED_BP} -lt "!{params.min_basepairs}" ]; then
            ERROR=1
            echo "!{sample} FASTQ(s) contain ${SEQUENCED_BP} total basepairs. This does not
                exceed the required minimum !{params.min_basepairs} bp. Further analysis is
                discontinued." | \
            sed 's/^\s*//' > !{sample}-low-sequence-depth-error.txt
        fi

        if [ ${TOTAL_READS} -lt "!{params.min_reads}" ]; then
            ERROR=1
            echo "!{sample} FASTQ(s) contain ${TOTAL_READS} total reads. This does not
                exceed the required minimum !{params.min_reads} read count. Further analysis is
                discontinued." | \
            sed 's/^\s*//' > !{sample}-low-read-count-error.txt
        fi
    fi

    if [ "${ERROR}" -eq "0" ]; then
        mkdir -p fastqs
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
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
fi
