// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'pirate')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process PIRATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::pirate=1.0.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate%3A1.0.4--hdfd78af_2' :
        'quay.io/biocontainers/pirate:1.0.4--hdfd78af_2' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                , emit: results
    tuple val(meta), path("core-genome.aln.gz")       , emit: aln
    tuple val(meta), path("gene_presence_absence.csv"), emit: csv
    path "*.{stdout.txt,stderr.txt,log,err}"          , emit: logs, optional: true
    path ".command.*"                                 , emit: nf_logs
    path "versions.yml"                               , emit: versions

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
    PIRATE_to_roary.pl -i results/ -o results/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}
    cp results/core_alignment.fasta.gz ./core-genome.aln.gz
    cp results/gene_presence_absence.csv ./
    touch gene_presence_absence.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
