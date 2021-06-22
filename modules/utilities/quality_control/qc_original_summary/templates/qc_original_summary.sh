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
        check-staging.py --fq1 !{fq[0]} --extra !{extra} --genome_size !{genome_size} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{extra} --genome_size !{genome_size}
    fi
fi

GENOME_SIZE=`head -n 1 !{genome_size}`
if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R1-original.json
    gzip -cd !{fq[1]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R2-original.json
    ln -s !{fq[0]} !{sample}_R1-original.fastq.gz
    ln -s !{fq[1]} !{sample}_R2-original.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}_R1-original.fastq.gz !{sample}_R2-original.fastq.gz
else
    # Single-End Reads
    gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}-original.json
    ln -s !{fq[0]} !{sample}-original.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}-original.fastq.gz
fi

mkdir -p quality-control/summary-original
mv *.json  quality-control/summary-original
mv *fastqc.html quality-control/summary-original
mv *fastqc.zip quality-control/summary-original

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
