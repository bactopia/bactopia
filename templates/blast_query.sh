#!/bin/bash
set -e
set -u

type=`readlink -f !{query}`
if [[ ${type} == *"blast/primers"* ]]; then
    mkdir -p blast/primers
    echo "#outfmt:!{params.outfmt}" > blast/primers/!{query_name}.txt
    blastn -db !{sample} \
           -query !{query} \
           -dust no \
           -word_size 7 \
           -perc_identity 100 \
           -evalue 1 \
           -num_threads !{task.cpus} \
           -outfmt '!{params.outfmt}' >> blast/primers/!{query_name}.txt
    #pigz --best -n -p !{task.cpus} blast/primers/!{query_name}.txt
elif [[ ${type} == *"blast/proteins"* ]]; then
    mkdir -p blast/proteins
    echo "#outfmt:!{params.outfmt}" > blast/proteins/!{query_name}.txt
    tblastn -db !{sample} \
            -query !{query} \
            -evalue 0.0001 \
            -num_threads !{task.cpus} \
            -outfmt '!{params.outfmt}' \
            -qcov_hsp_perc !{params.qcov_hsp_perc} >> blast/proteins/!{query_name}.txt
    #pigz --best -n -p !{task.cpus} blast/proteins/!{query_name}.txt
else
    mkdir -p blast/genes
    echo "#outfmt:!{params.outfmt}" > blast/genes/!{query_name}.txt
    blastn -db !{sample} \
           -query !{query} \
           -task blastn \
           -num_threads !{task.cpus} \
           -evalue 1 \
           -perc_identity !{params.perc_identity} \
           -qcov_hsp_perc !{params.qcov_hsp_perc} \
           -outfmt '!{params.outfmt}' >> blast/genes/!{query_name}.txt
    #pigz --best -n -p !{task.cpus} blast/genes/!{query_name}.txt
fi
