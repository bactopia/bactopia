// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'shigeifinder')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::shigeifinder=1.3.2"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SHIGEIFINDER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigeifinder:1.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/shigeifinder:1.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi
    shigeifinder \\
        $options.args \\
        --output ${prefix}.tsv \\
        -t $task.cpus \\
        -i $fasta_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: $VERSION
    END_VERSIONS
    """
}
