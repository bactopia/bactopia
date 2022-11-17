// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'mykrobe')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::mykrobe=0.12.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process MYKROBE_PREDICT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mykrobe:0.12.1--py39h70e0db4_1' :
        'quay.io/biocontainers/mykrobe:0.12.1--py39h70e0db4_1' }"

    input:
    tuple val(meta), path(seqs)
    val species

    output:
    tuple val(meta), path("${prefix}.csv") , emit: csv
    tuple val(meta), path("${prefix}.json"), emit: json
    path "*.{log,err}"     , optional: true, emit: logs
    path ".command.*"                      , emit: nf_logs
    path "versions.yml"                    , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    is_ont = meta.runtype == "ont" ? "--ont" : ""
    """
    mykrobe \\
        predict \\
        $options.args $is_ont \\
        --species $species \\
        --threads $task.cpus \\
        --sample $prefix \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
