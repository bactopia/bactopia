// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'tblastn')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::blast=2.14.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BLAST_TBLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.0--h7d5a4b4_1' :
        'quay.io/biocontainers/blast:2.14.0--h7d5a4b4_1' }"

    input:
    tuple val(meta), path(blastdb)
    path query

    output:
    tuple val(meta), path('*.tblastn.tsv'), emit: tsv
    path "*.{log,err}"                   , emit: logs, optional: true
    path ".command.*"                    , emit: nf_logs
    path "versions.yml"                  , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    // genes -> ffn, contigs -> fna
    def db_type = params.tblastn_use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${params.tblastn_outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf $blastdb
    
    ${which_cat} $query | \\
    tblastn \\
        -num_threads $task.cpus \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        $options.args \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "$outcols" | sed 's/<TAB>/\t/g' > ${prefix}.tblastn.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
