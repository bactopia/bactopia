#!/bin/bash
set -e
set -u

OUTDIR=primers
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/blast_primers.dry_run.txt
else
  for fasta in *.fasta; do
      type=`readlink -f ${fasta}`
      name="${fasta%.*}"
      mkdir -p ${OUTDIR}
      cat ${fasta} |
      parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
      blastn -db !{sample} \
             -outfmt \'!{params.outfmt}\' \
             -dust no \
             -word_size 7 \
             -perc_identity !{params.perc_identity} \
             -evalue 1 \
             -query - > ${OUTDIR}/${name}.txt

      if [[ !{params.compress} == "true" ]]; then
          pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.txt
      fi
  done
fi
