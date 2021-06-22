#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

echo "# FastQC Version" >> ${LOG_DIR}/!{task.process}.versions
fastqc -version>> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# fastq-scan Version" >> ${LOG_DIR}/!{task.process}.versions
fastq-scan -v >> ${LOG_DIR}/!{task.process}.versions 2>&1

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --genome_size !{genome_size} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --genome_size !{genome_size}
    fi
fi

GENOME_SIZE=`head -n 1 !{genome_size}`
if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R1-final.json
    gzip -cd !{fq[1]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R2-final.json
    ln -s !{fq[0]} !{sample}_R1-final.fastq.gz
    ln -s !{fq[1]} !{sample}_R2-final.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}_R1-final.fastq.gz !{sample}_R2-final.fastq.gz
else
    # Single-End Reads
    gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}-final.json
    ln -s !{fq[0]} !{sample}-final.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}-final.fastq.gz
fi

mkdir -p quality-control/summary-final
mv *.json  quality-control/summary-final
mv *fastqc.html quality-control/summary-final
mv *fastqc.zip quality-control/summary-final

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
