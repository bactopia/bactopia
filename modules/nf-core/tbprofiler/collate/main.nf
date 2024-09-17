// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'tbprofiler')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::tb-profiler=6.3.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:6.3.0--pyhdfd78af_0' :
        'quay.io/biocontainers/tb-profiler:6.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(json, stageAs: 'results/*')

    output:
    tuple val(meta), path("tbprofiler.csv")   , emit: csv
    tuple val(meta), path("tbprofiler.variants.csv"), emit: variants_csv
    tuple val(meta), path("tbprofiler.variants.txt"), emit: variants_txt
    tuple val(meta), path("*.itol..txt")            , emit: itol, optional: true
    path "*.{log,err}"                              , emit: logs, optional: true
    path ".command.*"                               , emit: nf_logs
    path "versions.yml"                             , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    tb-profiler collate --help

    tb-profiler \\
        collate \\
        $options.args \\
        --format csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler version 2>&1) | sed 's/tb-profiler version //')
    END_VERSIONS
    """
}
