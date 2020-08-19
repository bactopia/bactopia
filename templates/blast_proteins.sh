#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
OUTDIR=proteins
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/blast_proteins.dry_run.json
else
    mkdir -p ${LOG_DIR}
    touch ${LOG_DIR}/!{task.process}.versions
    echo "# tblastn Version" >> ${LOG_DIR}/!{task.process}.versions
    tblastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    echo "# Parallel Version" >> ${LOG_DIR}/!{task.process}.versions
    parallel --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    for fasta in *.fasta; do
        type=`readlink -f ${fasta}`
        name="${fasta%.*}"
        mkdir -p ${OUTDIR} temp_json
        cat ${fasta} |
        parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
        tblastn -db !{sample} \
                -outfmt 15 \
                -evalue 0.0001 \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
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
        cp .command.run ${LOG_DIR}/!{task.process}.run
        cp .command.sh ${LOG_DIR}/!{task.process}.sh
        cp .command.trace ${LOG_DIR}/!{task.process}.trace
    else
        rm -rf ${LOG_DIR}/
    fi
fi
