// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'pirate')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::pirate=1.0.4"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PIRATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate%3A1.0.4--hdfd78af_2' :
        'quay.io/biocontainers/pirate:1.0.4--hdfd78af_2' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("results/core_alignment.fasta.gz")  , emit: aln, optional: true
    tuple val(meta), path("results/gene_presence_absence.csv"), emit: csv, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    find . -name "*.gff.gz" | sed 's/.gz//' | xargs -I {} bash -c 'gzip -cdf {}.gz > {}'

    PIRATE \\
        $options.args \\
        --align \\
        --threads $task.cpus \\
        --input ./ \\
        --output results/
    PIRATE_to_roary.pl -i results/PIRATE.*.tsv -o results/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
