#!/bin/bash
set -e
set -u

mkdir blastdb
zcat !{fasta} | \
makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
