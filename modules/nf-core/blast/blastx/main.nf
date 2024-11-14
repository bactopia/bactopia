// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'blastx')
options.btype = "tools"
conda_tools   = "bioconda::blast=2.16.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BLAST_BLASTX {
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
    tuple val(meta), path('*.blastx.tsv'), emit: tsv
    path "*.{log,err}"                   , emit: logs, optional: true
    path ".command.*"                    , emit: nf_logs
    path "versions.yml"                  , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${params.blastx_outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf $blastdb
    
    ${which_cat} $query | \\
    blastx \\
        -num_threads $task.cpus \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.faa \\
        -query - \\
        $options.args \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "$outcols" | sed 's/<TAB>/\t/g' > ${prefix}.blastx.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.blastx.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastx: \$(blastx -version 2>&1 | sed 's/^.*blastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
