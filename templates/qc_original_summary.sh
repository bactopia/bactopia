#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
touch ${LOG_DIR}/!{task.process}.versions
echo "# FastQC Version" >> ${LOG_DIR}/!{task.process}.versions
fastqc -version>> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# fastq-scan Version" >> ${LOG_DIR}/!{task.process}.versions
fastq-scan -v >> ${LOG_DIR}/!{task.process}.versions 2>&1

GENOME_SIZE=`head -n 1 !{genome_size}`
if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    zcat !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R1-original.json
    zcat !{fq[1]} | fastq-scan -g ${GENOME_SIZE} > !{sample}_R2-original.json
    ln -s !{fq[0]} !{sample}_R1-original.fastq.gz
    ln -s !{fq[1]} !{sample}_R2-original.fastq.gz
    fastqc --noextract -f fastq -t !{task.cpus} !{sample}_R1-original.fastq.gz !{sample}_R2-original.fastq.gz
else
    # Single-End Reads
    zcat !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > !{sample}-original.json
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
    cp .command.run ${LOG_DIR}/!{task.process}.run
    cp .command.sh ${LOG_DIR}/!{task.process}.sh
    cp .command.trace ${LOG_DIR}/!{task.process}.trace
else
    rm -rf ${LOG_DIR}/
fi
