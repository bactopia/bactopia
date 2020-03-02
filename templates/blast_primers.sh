#!/bin/bash
set -e
set -u

for fasta in *.fasta; do
    type=`readlink -f ${fasta}`
    name="${fasta%.*}"
    mkdir -p primers
    cat ${fasta} |
    parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
    blastn -db !{sample} \
           -outfmt \'!{params.outfmt}\' \
           -dust no \
           -word_size 7 \
           -perc_identity !{params.perc_identity} \
           -evalue 1 \
           -query - > primers/${name}.txt

    if [[ !{params.compress} == "true" ]]; then
        pigz -n --best -p !{task.cpus} primers/${name}.txt
    fi
done
