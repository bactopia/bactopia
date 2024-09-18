// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'blastn')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::blast=2.16.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_2' :
        'quay.io/biocontainers/blast:2.16.0--hc155240_2' }"

    input:
    tuple val(meta), path(blastdb)
    path query

    output:
    tuple val(meta), path('*.blastn.tsv'), emit: tsv
    path "*.{log,err}"                   , emit: logs, optional: true
    path ".command.*"                    , emit: nf_logs
    path "versions.yml"                  , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    // genes -> ffn, contigs -> fna
    def db_type = params.blastn_use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${params.blastn_outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf $blastdb
    
    ${which_cat} $query | \\
    blastn \\
        -num_threads $task.cpus \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        $options.args \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "$outcols" | sed 's/<TAB>/\t/g' > ${prefix}.blastn.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.blastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
