// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'pbptyper')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::pbptyper=1.0.4" 
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PBPTYPER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}",  mode: params.publish_dir_mode,  overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbptyper:1.0.4--hdfd78af_0' :
        'quay.io/biocontainers/pbptyper:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("*.tblastn.tsv"), emit: blast
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbptyper \\
        $options.args \\
        --prefix $prefix \\
        --assembly $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/^.*pbptyper, version //;' )
    END_VERSIONS
    """
}
