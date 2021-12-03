// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'roary')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process ROARY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::roary=3.13.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0' :
        'quay.io/biocontainers/roary:3.13.0--pl526h516909a_0' }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("results/core-genome.aln.gz")      , emit: aln
    tuple val(meta), path("results/gene_presence_absence.csv"), emit: csv
    path "*.{stdout.txt,stderr.txt,log,err}"                  , emit: logs, optional: true
    path ".command.*"                                         , emit: nf_logs
    path "versions.yml"                                       , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mkdir gff
    cp -P gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs gunzip
    roary \\
        $options.args \\
        -p $task.cpus \\
        -f results/ \\
        gff/*.gff

    mv results/core_gene_alignment.aln results/core-genome.aln
    gzip results/*.aln
    gzip results/*.fa

    cat <<-END_VERSIONS > versions.yml
    roary:
        roary: \$( roary --version )
    END_VERSIONS
    """
}
