#!/bin/bash
set -e
set -u
OUTDIR=!{method}
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" >> ${LOG_DIR}/!{task.process}-!{method}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}-!{method}.versions

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

# Verify AWS files were staged
if [[ ! -L "!{fasta}" ]]; then
    check-staging.py --assembly !{fasta} --genome_size !{genome_size}
fi

if [ "!{method}" == "checkm" ]; then
    # CheckM
    mkdir checkm/
    if [ "$(uname)" = Darwin ]; then
        echo "checkm is not available due to pplacer not being available on MacOSX (via BioConda)" > checkm/checkm-not-available-on-macosx.txt
    elif [[ "!{params.skip_checkm}" == "true" ]]; then
        echo "checkm was skipped due to '--skip_checkm'" > checkm/checkm-was-skipped.txt
    else
        echo "# CheckM Version" >> ${LOG_DIR}/!{task.process}-!{method}.versions
        checkm -h | grep ":::" >> ${LOG_DIR}/!{task.process}-!{method}.versions 2>&1

        checkm lineage_wf ./ checkm/ \
            !{full_tree} --alignment_file checkm/checkm-genes.aln \
            --tab_table \
            --file checkm/checkm-results.txt \
            --threads !{task.cpus} \
            !{checkm_ali} !{checkm_nt}  --pplacer_threads !{task.cpus} \
            !{force_domain} !{no_refinement} --unique !{params.checkm_unique} \
            !{individual_markers} !{skip_adj_correction}  --multi !{params.checkm_multi} \
            !{skip_pseudogene_correction} !{ignore_thresholds} --aai_strain !{params.aai_strain} \
            --length !{params.checkm_length} > ${LOG_DIR}/checkm.out 2> ${LOG_DIR}/checkm.err

        if [[ !{params.compress} == "true" ]]; then
            find . -name "*.faa" -or -name "*hmmer.analyze.txt" | xargs -I {} pigz -n --best -p !{task.cpus} {}
        fi
    fi
else
    # QUAST
    echo "# QUAST Version" >> ${LOG_DIR}/!{task.process}-!{method}.versions
    quast --version >> ${LOG_DIR}/!{task.process}-!{method}.versions 2>&1
    GENOME_SIZE=`head -n 1 !{genome_size}`
    est_ref_size=""
    if [ "${GENOME_SIZE}" != "0" ]; then
        est_ref_size="--est-ref-size ${GENOME_SIZE}"
    fi
    quast !{fasta} ${est_ref_size} \
        -o quast \
        --threads !{task.cpus} \
        --glimmer \
        --contig-thresholds !{params.contig_thresholds} \
        --plots-format !{params.plots_format} > ${LOG_DIR}/quast.out 2> ${LOG_DIR}/quast.err
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}-!{method}.err
    cp .command.out ${LOG_DIR}/!{task.process}-!{method}.out
    cp .command.sh ${LOG_DIR}/!{task.process}-!{method}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}-!{method}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
