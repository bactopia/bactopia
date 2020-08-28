#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
if [ "!{params.dry_run}" == "true" ]; then
    touch !{sample}.sig !{sample}-k31.msh !{sample}-k21.msh
else
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
    mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
    sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    gzip -cd !{fastq} | mash sketch -o !{sample}-k21 -k 21 -s !{params.mash_sketch} -r -I !{sample} -
    gzip -cd !{fastq} | mash sketch -o !{sample}-k31 -k 31 -s !{params.mash_sketch} -r -I !{sample} -
    sourmash compute --scaled !{params.sourmash_scale} -o !{sample}.sig -p !{task.cpus} \
                     --track-abundance --merge !{sample} -k 21,31,51 !{fastq}

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
fi
