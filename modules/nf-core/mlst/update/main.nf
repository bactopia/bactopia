// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options     = initOptions(params.containsKey("options") ? params.options : [:], 'mlst_update')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::mlst=2.23.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process MLST_UPDATE {
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_1' :
        'quay.io/biocontainers/mlst:2.23.0--hdfd78af_1' }"

    output:
    path "mlst.tar.gz" , emit: db
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = "mlst-update"
    """
    # Download Data
    mkdir -p mlstdb/pubmlst mlstdb/blast
    mlst-download_pub_mlst -j ${task.cpus} -d mlstdb/pubmlst

    # Make BLAST db
    # Modified from https://github.com/tseemann/mlst/blob/master/scripts/mlst-make_blast_db
    for N in \$(find mlstdb/pubmlst -mindepth 1 -maxdepth 1 -type d); do
        SCHEME=\$(basename \$N)
        echo "Adding: \$SCHEME"
        cat mlstdb/pubmlst/\$SCHEME/*.tfa | grep -v 'not a locus' | sed -e "s/^>/>\$SCHEME./" >> mlstdb/blast/mlst.fa
    done
    makeblastdb -hash_index -in mlstdb/blast/mlst.fa -dbtype nucl -title "PubMLST" -parse_seqids
    echo "Created BLAST database for mlstdb/blast/mlst.fa"

    # Create MLST database version file
    MLST_VERSION=\$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    MLSTDB_VERSION=\$(date +"%Y%m%d")
    echo "\${MLST_VERSION}-\${MLSTDB_VERSION}" > mlstdb/DB_VERSION

    # Wrap it all up
    tar -czvf mlst.tar.gz mlstdb/ 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
        mlst-database: \$( cat mlstdb/DB_VERSION )
    END_VERSIONS
    """
}
