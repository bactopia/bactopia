#!/bin/bash
set -e
set -u

OUTDIR=proteins
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/blast_proteins.dry_run.txt
else
    for fasta in *.fasta; do
        type=`readlink -f ${fasta}`
        name="${fasta%.*}"
        mkdir -p ${OUTDIR}
        cat ${fasta} |
        parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
        tblastn -db !{sample} \
                -outfmt \'!{params.outfmt}\' \
                -evalue 0.0001 \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - > ${OUTDIR}/${name}.txt

        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.txt
        fi
    done
fi
