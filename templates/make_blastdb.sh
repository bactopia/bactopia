#!/bin/bash
set -e
set -u

mkdir blastdb
if [[ !{params.compress} == "true" ]]; then
    zcat !{fasta} | \
    makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
else
    cat !{fasta} | \
    makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
fi
