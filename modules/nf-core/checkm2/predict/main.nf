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
    tuple val(meta), path("results/*"),                             emit: results
    tuple val(meta), path("results/quality_report.tsv"),  emit: tsv
    path "*.{log,err}",                                             emit: logs, optional: true
    path ".command.*",                                              emit: nf_logs
    path "versions.yml",                                            emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def db_is_compressed = db.getName().endsWith(".dmnd.gz") ? true : false
    def db_is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    echo $db
    
    # Decompress fasta file if compressed
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Decompress fasta file if compressed
    if [ "$db_is_compressed" == "true" ]; then
        mkdir database
        gunzip -c $db > database/checkm2_db.dmnd
        CHECKM2_DB=\$(find database/ -name "*.dmnd")
        fi

    # Check if db is a tarball and if so extract and decompress it
    if [ "$db_is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        CHECKM2_DB=\$(find database/ -name "*.dmnd")
    fi

    # Check if db is a directory
    if [ -d "$db" ]; then
        echo "Database is a directory, expected a .dmnd file. Searching for .dmnd file in directory"
        CHECKM2_DB=\$(find ${db}/ -name "*.dmnd")
    # Else Check if db ends in .dmnd
    elif [ "$db" == "*.dmnd" ]; then
        CHECKM2_DB=$db
    fi

    # Check if CHECKM2_DB is set
    if [ -z "\$CHECKM2_DB" ]; then
        echo "ERROR: No database found. Please provide a .dmnd file or a directory containing a .dmnd file"
        exit 1
    fi

    checkm2 \\
        predict \\
        --input ${fasta} \\
        --output-directory results \\
        --threads ${task.cpus} \\
        --database_path \$CHECKM2_DB \\
        $options.args

    mv results/checkm2.log ./

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
