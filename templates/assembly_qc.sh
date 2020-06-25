#!/bin/bash
set -e
set -u
OUTDIR=!{method}
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/!{method}.dry_run.txt
else

    GENOME_SIZE=`head -n 1 !{genome_size}`
    if [ "!{method}" == "checkm" ]; then
        # CheckM
    else
        # QUAST
    fi

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate files
    fi
fi
