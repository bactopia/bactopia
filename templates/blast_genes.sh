#!/bin/bash
set -e
set -u

OUTDIR=genes
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/blast_genes.dry_run.txt
else
    for fasta in *.fasta; do
        type=`readlink -f ${fasta}`
        name="${fasta%.*}"
        mkdir -p ${OUTDIR}
        cat ${fasta} |
        parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
        blastn -db !{sample} \
            -outfmt \'!{params.outfmt}\' \
            -task blastn \
            -evalue 1 \
            -perc_identity !{params.perc_identity} \
            -qcov_hsp_perc !{params.qcov_hsp_perc} \
            -query - > ${OUTDIR}/${name}.txt

        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.txt
        fi
    done
fi
