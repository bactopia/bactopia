#!/bin/bash
set -e
set -u

if [ "!{method}" == "blast" ]; then
    mkdir blast
    !{baseDir}/bin/mlst-blast.py !{assembly} !{database} blast/!{sample}-blast.json \
        --cpu !{task.cpus} --compressed
elif [ "!{method}" == "ariba" ]; then
    ariba run !{database} !{fq[0]} !{fq[1]} ariba --threads !{task.cpus} --verbose
fi
