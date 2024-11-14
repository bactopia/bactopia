// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'tblastx')
options.btype = "tools"
conda_tools   = "bioconda::blast=2.16.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BLAST_TBLASTX {
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
    tuple val(meta), path('*.tblastx.tsv'), emit: tsv
    path "*.{log,err}"                   , emit: logs, optional: true
    path ".command.*"                    , emit: nf_logs
    path "versions.yml"                  , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    // genes -> ffn, contigs -> fna
    def db_type = params.tblastx_use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${params.tblastx_outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf $blastdb
    
    ${which_cat} $query | \\
    tblastx \\
        -num_threads $task.cpus \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        $options.args \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "$outcols" | sed 's/<TAB>/\t/g' > ${prefix}.tblastx.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastx.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastx: \$(tblastx -version 2>&1 | sed 's/^.*tblastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
