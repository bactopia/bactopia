#!/bin/bash
set -e
set -u

zcat !{genes} | \
blastn -db !{blastdb} -out !{sample}-plsdb.json -outfmt '15' -task blastn \
       -evalue 1 -num_threads !{task.cpus}
pigz --best -n -p !{task.cpus} !{sample}-plsdb.json
