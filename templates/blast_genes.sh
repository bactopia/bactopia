#!/bin/bash
set -e
set -u

for fasta in *.fasta; do
    type=`readlink -f ${fasta}`
    name="${fasta%.*}"
    mkdir -p genes
    cat ${fasta} |
    parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
    blastn -db !{sample} \
        -outfmt \'!{params.outfmt}\' \
        -task blastn \
        -evalue 1 \
        -perc_identity !{params.perc_identity} \
        -qcov_hsp_perc !{params.qcov_hsp_perc} \
        -query - > genes/${name}.txt

    if [[ !{params.compress} == "true" ]]; then
        pigz -n --best -p !{task.cpus} genes/${name}.txt
    fi
done
