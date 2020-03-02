#!/bin/bash
set -e
set -u

for fasta in *.fasta; do
    type=`readlink -f ${fasta}`
    name="${fasta%.*}"
    mkdir -p proteins
    cat ${fasta} |
    parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
    tblastn -db !{sample} \
            -outfmt \'!{params.outfmt}\' \
            -evalue 0.0001 \
            -qcov_hsp_perc !{params.qcov_hsp_perc} \
            -query - > proteins/${name}.txt

    if [[ !{params.compress} == "true" ]]; then
        pigz -n --best -p !{task.cpus} proteins/${name}.txt
    fi
done
