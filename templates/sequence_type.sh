#!/bin/bash
set -e
set -u

tar -xzvf !{dataset_tarball}

if [ "!{method}" == "blast" ]; then
    mkdir -p blast
    mlst-blast.py !{assembly} !{dataset_name} blast/!{sample}-blast.json \
        --cpu !{task.cpus} --compressed
elif [ "!{method}" == "ariba" ]; then
    mv !{dataset_name}/ref_db ./
    ariba run ref_db !{fq[0]} !{fq[1]} ariba \
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
