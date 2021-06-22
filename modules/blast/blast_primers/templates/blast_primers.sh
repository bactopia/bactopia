#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
OUTDIR=primers
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# blastn Version" >> ${LOG_DIR}/!{task.process}.versions
blastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# Parallel Version" >> ${LOG_DIR}/!{task.process}.versions
parallel --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
for fasta in *.fasta; do
    type=`readlink -f ${fasta}`
    name="${fasta%.*}"
    mkdir -p ${OUTDIR} temp_json
    cat ${fasta} | sed -e 's/<[^>]*>//g' |
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

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
