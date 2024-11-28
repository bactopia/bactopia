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
    tuple val(meta), path(fasta)    //, stageAs: "input_bins/*")  // left from original nf-core module
    path db
    //tuple val(dbmeta), path(db) // left from original nf-core module

    output:
    tuple val(meta), path("${prefix}"),                             emit: results
    tuple val(meta), path("${prefix}/quality_report.tsv"),  emit: tsv
    path "*.{log,err}",                                             emit: logs, optional: true
    path ".command.*",                                              emit: nf_logs
    path "versions.yml",                                            emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    #if [ "$is_tarball" == "true" ]; then
    #    mkdir database
    #    tar -xzf $db -C database
    #    CHECKM2_DB=\$(find database/ -name "*.dmnd")
    #else
    #    CHECKM2_DB=\$db
    #fi

    checkm2 \\
        predict \\
        --input ${fasta} \\
        --output-directory ${prefix} \\
        --threads ${task.cpus} \\
        --database_path $db \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/diamond_output ${prefix}/protein_files
    touch ${prefix}/quality_report.tsv ${prefix}/checkm2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
