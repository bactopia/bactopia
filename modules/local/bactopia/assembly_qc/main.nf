nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'assembly_qc')

process ASSEMBLY_QC {
    /* Assess the quality of the assembly using QUAST and CheckM */
    tag "${meta.id}"
    label "max_cpu_75"
    label "process_medium"
    label "assembly_qc"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(genome_size), path(fasta), path(total_contigs)

    output:
    path "results/*"
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

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
    is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    '''
    if [ "!{is_compressed}" == "true" ]; then
        gzip -c -d !{fasta} > !{fasta_name}
    fi

    RAN_CHECKM="0"
    # CheckM
    mkdir checkm/
    if [ "$(uname)" = Darwin ]; then
        echo "checkm is not available due to pplacer not being available on MacOSX (via BioConda)" > checkm/checkm-not-available-on-macosx.txt
    elif [[ "!{params.run_checkm}" == "true" ]]; then
        RAN_CHECKM="1"
        checkm lineage_wf ./ checkm/ \
            --alignment_file checkm/checkm-genes.aln \
            --tab_table \
            --file checkm/checkm-results.txt \
            --threads !{task.cpus} \
            --pplacer_threads !{task.cpus} \
            --unique !{params.checkm_unique} \
            --multi !{params.checkm_multi} \
            --aai_strain !{params.aai_strain} \
            --length !{params.checkm_length} !{checkm_opts} > checkm.stdout.txt 2> checkm.stderr.txt
        mv checkm/checkm.log ./

        if [[ !{params.skip_compression} == "false" ]]; then
            find . -name "*.faa" -or -name "*hmmer.analyze.txt" | xargs -I {} pigz -n --best -p !{task.cpus} {}
        fi
    else
        echo "checkm was skipped due to not including '--run_checkm'" > checkm/checkm-was-skipped.txt
    fi

    # QUAST
    GENOME_SIZE=`head -n 1 !{genome_size}`
    est_ref_size=""
    if [ "${GENOME_SIZE}" != "0" ]; then
        est_ref_size="--est-ref-size ${GENOME_SIZE}"
    fi

    quast !{fasta_name} ${est_ref_size} \
        -o quast \
        --threads !{task.cpus} \
        --glimmer \
        --contig-thresholds !{params.contig_thresholds} \
        --plots-format !{params.plots_format} > quast.stdout.txt 2> quast.stderr.txt
    mv quast/quast.log ./

    # Results dir
    mkdir results
    mv checkm/ results/
    mv quast/ results/

    # Capture versions (no indent or it'll break)
    if [ "${RAN_CHECKM}" == "1" ]; then

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        checkm: $(echo $(checkm -h 2>&1) | sed 's/.*CheckM v//;s/ .*$//')
        quast: $(echo $(quast --version 2>&1) | sed 's/.*QUAST v//;s/ .*$//')
    END_VERSIONS

    else

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        quast: $(echo $(quast --version 2>&1) | sed 's/.*QUAST v//;s/ .*$//')
    END_VERSIONS

    fi
    '''
}
