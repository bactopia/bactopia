#!/bin/bash
set -e
set -u

if [ "!{minmer_database}" == "refseq-k21-s1000.msh" ]; then
    zcat !{fastq} | \
    mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{database}  - | \
    sort -gr > !{sample}-refseq-k21.txt
elif [ "!{minmer_database}" == "genbank-k21.json.gz" ]; then
    sourmash lca gather !{sourmash} !{database} > !{sample}-genbank-k21.txt
elif [ "!{minmer_database}" == "genbank-k31.json.gz" ]; then
    sourmash lca gather !{sourmash} !{database} > !{sample}-genbank-k31.txt
else
    sourmash lca gather !{sourmash} !{database} > !{sample}-genbank-k51.txt
fi
