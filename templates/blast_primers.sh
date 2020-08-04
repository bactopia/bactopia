#!/bin/bash
set -e
set -u

OUTDIR=primers
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/blast_primers.dry_run.json
else
  for fasta in *.fasta; do
      type=`readlink -f ${fasta}`
      name="${fasta%.*}"
      mkdir -p ${OUTDIR} temp_json
      cat ${fasta} |
      parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
      blastn -db !{sample} \
             -outfmt 15 \
             -task blastn \
             -dust no \
             -word_size 7 \
             -perc_identity !{params.perc_identity} \
             -evalue 1 \
             -query - \
             -out temp_json/${name}_{#}.json

      merge-blast-json.py temp_json > ${OUTDIR}/${name}.json
      rm -rf temp_json

      if [[ !{params.compress} == "true" ]]; then
          pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.json
      fi
  done
fi
