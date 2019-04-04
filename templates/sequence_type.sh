#!/bin/bash
set -e
set -u

if [ "!{method}" == "blast" ]; then
    mkdir blast
    !{baseDir}/bin/mlst-blast.py !{assembly} !{dataset} blast/!{sample}-blast.json \
        --cpu !{task.cpus} --compressed
elif [ "!{method}" == "ariba" ]; then
    ariba run !{dataset} !{fq[0]} !{fq[1]} ariba \
        --nucmer_min_id !{params.nucmer_min_id} \
        --nucmer_min_len !{params.nucmer_min_len} \
        --nucmer_breaklen !{params.nucmer_breaklen} \
        --assembly_cov !{params.assembly_cov} \
        --min_scaff_depth !{params.min_scaff_depth} \
        --assembled_threshold !{params.assembled_threshold} \
        --gene_nt_extend !{params.gene_nt_extend} \
        --unique_threshold !{params.unique_threshold} \
        --threads !{task.cpus} \
        --force \
        --noclean \
        --verbose !{spades_options}
fi
