#!/bin/bash
set -e
set -u

file_size=`gzip -dc !{genes} | wc -c`
block_size=$(( file_size / !{task.cpus} / 2 ))
zcat !{genes} | \
parallel --gnu --plain -j !{task.cpus} --block ${block_size} --recstart '>' --pipe \
blastn -db !{blastdb} \
       -outfmt 15 \
       -task blastn \
       -evalue 1 \
       -max_target_seqs !{params.max_target_seqs} \
       -perc_identity !{params.perc_identity} \
       -qcov_hsp_perc !{params.qcov_hsp_perc} \
       -query - > !{sample}-plsdb.json
