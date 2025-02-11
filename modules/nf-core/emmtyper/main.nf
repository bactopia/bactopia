// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'emmtyper')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::emmtyper=0.2.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process EMMTYPER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emmtyper:0.2.0--py_0' :
        'quay.io/biocontainers/emmtyper:0.2.0--py_0' }"

    input:
    tuple val(meta), path(fasta)
    path blastdb

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    echo $db

    # Conditionally add the database if it is provided by user
    if [ "$blastdb" == "" ]; then
        emmtyper \\
            $options.args \\
            $fasta_name \\
            > ${prefix}.tsv
    else

        # Make the blast database
        makeblastdb -in $blastdb -dbtype nucl

        emmtyper \\
            --blast_db $blastdb \\
            $options.args \\
            $fasta_name \\
            > ${prefix}.tsv
    fi

    # If 'tmp' is not in $fasta_name, remove '.tmp' from the output files contents
    if [ $fasta_name != *.tmp* ]; then
        sed -i 's/.tmp\t/\t/g' ${prefix}.tsv
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emmtyper: \$( echo \$(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' )
    END_VERSIONS
    """
}
