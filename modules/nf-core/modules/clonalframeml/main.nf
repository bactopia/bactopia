// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process CLONALFRAMEML {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::clonalframeml=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/clonalframeml:1.12--h7d875b9_1"
    } else {
        container "quay.io/biocontainers/clonalframeml:1.12--h7d875b9_1"
    }

    input:
    tuple val(meta), path(newick), path(msa)

    output:
    tuple val(meta), path("*.emsim.txt")                   , emit: emsim, optional: true
    tuple val(meta), path("*.em.txt")                      , emit: em
    tuple val(meta), path("*.importation_status.txt")      , emit: status
    tuple val(meta), path("*.labelled_tree.newick")        , emit: newick
    tuple val(meta), path("*.ML_sequence.fasta")           , emit: fasta
    tuple val(meta), path("*.position_cross_reference.txt"), emit: pos_ref
    path "*.{stdout.txt,stderr.txt,log,err}"               , emit: logs, optional: true
    path ".command.*"                                      , emit: nf_logs
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    ClonalFrameML \\
        $newick \\
        $msa_name \\
        $prefix \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
    END_VERSIONS
    """
}
