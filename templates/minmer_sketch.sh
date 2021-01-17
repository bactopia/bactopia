#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]}
    fi
fi

gzip -cd !{fastq} | mash sketch -o !{sample}-k21 -k 21 -s !{params.mash_sketch} -r -I !{sample} -
gzip -cd !{fastq} | mash sketch -o !{sample}-k31 -k 31 -s !{params.mash_sketch} -r -I !{sample} -
sourmash compute --scaled !{params.sourmash_scale} -o !{sample}.sig -p !{task.cpus} \
                    --track-abundance --merge !{sample} -k 21,31,51 !{fastq}

# pass the FASTQs along
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

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
