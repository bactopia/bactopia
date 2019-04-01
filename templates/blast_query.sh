#!/bin/bash
set -e
set -u

type=`readlink -f !{query}`
if [[ ${type} == *"blast/primers"* ]]; then
    mkdir -p blast/primers
    blastn -db !{sample} -query !{query} -dust no -word_size 7 -perc_identity 100 \
           -num_threads !{task.cpus} -outfmt 15 > blast/primers/!{query_name}.json
    pigz --best -n -p !{task.cpus} blast/primers/!{query_name}.json
elif [[ ${type} == *"blast/proteins"* ]]; then
    mkdir -p blast/proteins
    tblastn -db !{sample} -query !{query}  -evalue 0.0001 \
            -num_threads !{task.cpus} -outfmt 15 > blast/proteins/!{query_name}.json
    pigz --best -n -p !{task.cpus} blast/proteins/!{query_name}.json
else
    mkdir -p blast/genes
    blastn -db !{sample} -query !{query} -num_threads !{task.cpus} -outfmt 15 > blast/genes/!{query_name}.json
    pigz --best -n -p !{task.cpus} blast/genes/!{query_name}.json
fi
