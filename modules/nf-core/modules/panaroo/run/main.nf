// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'panaroo')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::panaroo=1.2.9"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PANAROO_RUN {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panaroo:1.2.9--pyhdfd78af_0' :
        'quay.io/biocontainers/panaroo:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                                              , emit: results
    tuple val(meta), path("core-genome.aln.gz")                     , optional: true, emit: aln
    tuple val(meta), path("results/gene_presence_absence_roary.csv"), optional: true, emit: csv
    tuple val(meta), path("results/gene_presence_absence.csv")      , optional: true, emit: panaroo_csv
    path "*.{log,err}"                                              , optional: true, emit: logs
    path ".command.*"                                                               , emit: nf_logs
    path "versions.yml"                                                             , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs gunzip

    panaroo \\
        $options.args \\
        -t $task.cpus \\
        -o results \\
        -i gff/*.gff

    # Cleanup
    find . -name "*.fas" | xargs -I {} -P $task.cpus -n 1 gzip {}

    if [[ -f "results/core_gene_alignment.aln" ]]; then
        gzip results/core_gene_alignment.aln
        cp results/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' ))
    END_VERSIONS
    """
}
