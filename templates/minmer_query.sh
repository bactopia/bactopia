#!/bin/bash
set -e
set -u

if [ "!{dataset_name}" == "refseq-k21-s1000.msh" ]; then
    zcat !{fastq} | \
    mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
    sort -gr > !{sample}-refseq-k21.txt
elif [ "!{dataset_name}" == "plsdb.msh" ]; then
    zcat !{fastq} | \
    mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
    sort -gr > !{sample}-plsdb-k21.txt
elif [ "!{dataset_name}" == "genbank-k21.json.gz" ]; then
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k21.txt
elif [ "!{dataset_name}" == "genbank-k31.json.gz" ]; then
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k31.txt
else
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k51.txt
fi
