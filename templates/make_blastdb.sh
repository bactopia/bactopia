#!/bin/bash
set -e
set -u

mkdir blastdb
if [ "!{params.dry_run}" == "true" ]; then
    touch blastdb/!{sample}.nin
else
    if [[ !{params.compress} == "true" ]]; then
        zcat !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample} -parse_seqids
    else
        cat !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample} -parse_seqids
    fi
fi
