// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process PIRATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::pirate=1.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pirate%3A1.0.4--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/pirate:1.0.4--hdfd78af_1"
    }

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                , emit: results
    tuple val(meta), path("core_alignment.fasta")     , emit: aln
    tuple val(meta), path("gene_presence_absence.csv"), emit: csv
    path "*.{stdout.txt,stderr.txt,log,err}"          , emit: logs, optional: true
    path ".command.*"                                 , emit: nf_logs
    path "versions.yml"                               , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    PIRATE \\
        $options.args \\
        --align \\
        --rplots \\
        --threads $task.cpus \\
        --input ./ \\
        --output results/
    PIRATE_to_roary.pl -i results.*.tsv -o results/gene_presence_absence.csv

    cp results/core_alignment.fasta ./
    gzip results/*.fasta
    cp results/gene_presence_absence.csv ./


    cat <<-END_VERSIONS > versions.yml
    pirate:
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
