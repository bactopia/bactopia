nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "assembly_qc"

process ASSEMBLY_QC {
    /* Assess the quality of the assembly using QUAST and CheckM */
    tag "${sample} - ${method}"
    label "max_cpu_75"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), path(genome_size), path(fasta), path(total_contigs)
    each method

    output:
    path "${method}/*"
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    //CheckM Related
    full_tree = params.full_tree ? '' : '--reduced_tree'
    checkm_ali = params.checkm_ali ? '--ali' : ''
    checkm_nt = params.checkm_nt ? '--nt' : ''
    force_domain = params.force_domain ? '--force_domain' : ''
    no_refinement = params.no_refinement ? '--no_refinement' : ''
    individual_markers = params.individual_markers ? '--individual_markers' : ''
    skip_adj_correction = params.skip_adj_correction ? '--skip_adj_correction' : ''
    skip_pseudogene_correction = params.skip_pseudogene_correction ? '--skip_pseudogene_correction' : ''
    ignore_thresholds = params.ignore_thresholds ? '--ignore_thresholds' : ''
    checkm_opts = [full_tree, checkm_ali, checkm_nt, force_domain, no_refinement, individual_markers, 
                   skip_adj_correction, skip_pseudogene_correction, ignore_thresholds].join(" ")
    '''
    if [ "!{method}" == "checkm" ]; then
        # CheckM
        mkdir checkm/
        if [ "$(uname)" = Darwin ]; then
            echo "checkm is not available due to pplacer not being available on MacOSX (via BioConda)" > checkm/checkm-not-available-on-macosx.txt
        elif [[ "!{params.skip_checkm}" == "true" ]]; then
            echo "checkm was skipped due to '--skip_checkm'" > checkm/checkm-was-skipped.txt
        else
            checkm lineage_wf ./ checkm/ \
                --alignment_file checkm/checkm-genes.aln \
                --tab_table \
                --file checkm/checkm-results.txt \
                --threads !{task.cpus} \
                --pplacer_threads !{task.cpus} \
                --unique !{params.checkm_unique} \
                --multi !{params.checkm_multi} \
                --aai_strain !{params.aai_strain} \
                --length !{params.checkm_length} !{checkm_opts} > checkm.stdout.txt 2> checkm.stderr.out

            if [[ !{params.compress} == "true" ]]; then
                find . -name "*.faa" -or -name "*hmmer.analyze.txt" | xargs -I {} pigz -n --best -p !{task.cpus} {}
            fi

            # Capture Version
            checkm -h | grep ":::" > checkm.version.txt 2>&1
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
            --plots-format !{params.plots_format} > quast.stdout.txt 2> quast.stderr.txt

        # Capture version
        quast --version > quast.version.txt 2>&1
    fi
    '''

    stub:
    """
    mkdir ${method}
    touch ${method}/${sample}
    """
}
