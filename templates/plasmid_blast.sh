#!/bin/bash
set -e
set -u

echo "#outfmt:!{params.outfmt}" > !{sample}-plsdb.txt
zcat !{genes} | \
blastn -db !{blastdb} \
       -outfmt '!{params.outfmt}' \
       -task blastn \
       -evalue 1 \
       -num_threads !{task.cpus} \
       -max_target_seqs !{params.max_target_seqs} \
       -perc_identity !{params.perc_identity} \
       -qcov_hsp_perc !{params.qcov_hsp_perc} >> !{sample}-plsdb.txt \
#pigz --best -n -p !{task.cpus} !{sample}-plsdb.txt
