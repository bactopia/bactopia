// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'checkm')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::checkm-genome=1.2.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.2.3--pyhdfd78af_0' :
        'quay.io/biocontainers/checkm-genome:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/*")    , emit: results
    tuple val(meta), path("results/${prefix}-results.txt"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    checkm \\
        lineage_wf ./ results/ \\
        --tab_table \\
        --threads $task.cpus \\
        --pplacer_threads $task.cpus \\
        --alignment_file results/${prefix}-genes.aln \\
        --file results/${prefix}-results.txt \\
        $options.args

    find ./results/ -name "*.faa" -or -name "*hmmer.analyze.txt" -or -name "*.fasta" | xargs gzip
    mv results/checkm.log ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(echo \$(checkm -h 2>&1) | sed 's/.*CheckM v//;s/ .*\$//')
    END_VERSIONS
    """
}
