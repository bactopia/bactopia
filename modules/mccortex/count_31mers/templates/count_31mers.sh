#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# mccortex31 Version" >> ${LOG_DIR}/!{task.process}.versions
mccortex31 2>&1 | grep "version" >> ${LOG_DIR}/!{task.process}.versions 2>&1

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]}
    fi
fi

if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    mccortex31 build -f -k 31 -s !{sample} -2 !{fq[0]}:!{fq[1]} -t !{task.cpus} -m !{m}mb -q temp_counts
else
    # Single-End Reads
    mccortex31 build -f -k 31 -s !{sample} -1 !{fq[0]} -t !{task.cpus} -m !{m}mb -q temp_counts
fi

if [ "!{params.keep_singletons}" == "false" ]; then
    # Clean up Cortex file (mostly remove singletons)
    mccortex31 clean -q -B 2 -U2 -T2 -m !{m}mb -o !{sample}.ctx temp_counts
    rm temp_counts
else
    mv temp_counts !{sample}.ctx
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
