// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'roary')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::roary=3.13.0" 
conda_env   = file("${params.condadir}/roary").exists() ? "${params.condadir}/roary" : conda_tools

process ROARY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0' :
        'quay.io/biocontainers/roary:3.13.0--pl526h516909a_0' }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("results/core-genome.aln.gz")       , emit: aln, optional: true
    tuple val(meta), path("results/gene_presence_absence.csv"), emit: csv, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs gunzip
    roary \\
        $options.args \\
        -p $task.cpus \\
        -f results/ \\
        gff/*.gff


    # Only copy files if they exist
    if [[ -f "results/core_alignment.fasta.gz" ]]; then
        cp results/core_alignment.fasta.gz ./core-genome.aln.gz
    fi

    cp results/core_gene_alignment.aln results/core-genome.aln
    gzip results/*.aln
    gzip results/*.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """
}
