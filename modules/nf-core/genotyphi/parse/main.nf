// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'genotyphi')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::genotyphi=1.9.1" 
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process GENOTYPHI_PARSE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genotyphi:1.9.1--hdfd78af_1' :
        'quay.io/biocontainers/genotyphi:1.9.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.tsv")     , emit: tsv
    path "*.{log,err}", optional: true , emit: logs
    path ".command.*"                  , emit: nf_logs
    path "versions.yml"                , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    parse_typhi_mykrobe.py \\
        --jsons $json \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genotyphi: \$(echo \$(genotyphi --version 2>&1) | sed 's/^.*GenoTyphi v//;' )
    END_VERSIONS
    """
}
