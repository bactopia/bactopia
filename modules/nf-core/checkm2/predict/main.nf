// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'checkm2')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::checkm2=1.0.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CHECKM2_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.2--pyh7cba7a3_0':
        'quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("results/*"),             emit: results
    tuple val(meta), path("results/${prefix}.tsv"), emit: tsv
    path "*.{log,err}",                             emit: logs, optional: true
    path ".command.*",                              emit: nf_logs
    path "versions.yml",                            emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    echo $db
    
    # Decompress fasta file if compressed
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Check if db is a directory - if so, find the diamond database
    if [ -d "$db" ]; then
        CHECKM2_DB=\$(find ${db}/ -name "*.dmnd")
    else
        CHECKM2_DB=$db
    fi

    checkm2 \\
        predict \\
        --input ${fasta} \\
        --output-directory results \\
        --threads ${task.cpus} \\
        --database_path \$CHECKM2_DB \\
        $options.args

    mv results/checkm2.log ./
    mv results/quality_report.tsv results/${prefix}.tsv

    # Cleanup
    gzip results/protein_files/*.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
