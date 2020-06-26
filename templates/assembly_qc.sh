#!/bin/bash
set -e
set -u
OUTDIR=!{method}
if [ "!{params.dry_run}" == "true" ]; then
    mkdir ${OUTDIR}
    touch ${OUTDIR}/!{method}.dry_run.txt
else
    if [ "!{method}" == "checkm" ]; then
        # CheckM
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
            --length !{params.checkm_length}

        if [[ !{params.compress} == "true" ]]; then
            find -name "*.faa" -or -name "*hmmer.analyze.txt" | xargs -I {} pigz -n --best -p !{task.cpus} {}
        fi
    else
        # QUAST
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
            --plots-format !{params.plots_format}
    fi
fi
