#!/bin/bash
set -e
set -u

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
