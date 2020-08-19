#!/bin/bash
set -e
set -u
OUTDIR=!{method}
LOG_DIR="!{task.process}"
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/!{method}.dry_run.txt
else
    mkdir -p ${LOG_DIR}
    if [ "!{method}" == "checkm" ]; then
        # CheckM
        echo "# CheckM Version" > ${LOG_DIR}/!{task.process}-!{method}.versions
        checkm -h | grep ":::" >> ${LOG_DIR}/!{task.process}-!{method}.versions 2>&1

        mkdir checkm/
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
            find -name "*.faa" -or -name "*hmmer.analyze.txt" | xargs -I {} pigz -n --best -p !{task.cpus} {}
        fi
    else
        # QUAST
        echo "# QUAST Version" > ${LOG_DIR}/!{task.process}-!{method}.versions
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
        cp .command.run ${LOG_DIR}/!{task.process}-!{method}.run
        cp .command.sh ${LOG_DIR}/!{task.process}-!{method}.sh
        cp .command.trace ${LOG_DIR}/!{task.process}-!{method}.trace
    else
        rm -rf ${LOG_DIR}/
    fi
fi
