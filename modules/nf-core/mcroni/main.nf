def VERSION = '1.0.4' // Version information not provided by tool on CLI
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'mcroni')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::mcroni=1.0.4 conda-forge::python=3.9"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process MCRONI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mcroni:1.0.4--pyh5e36f6f_0' :
        'quay.io/biocontainers/mcroni:1.0.4--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")               , emit: tsv
    tuple val(meta), path("*.fa"), optional: true, emit: fa
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mcroni \\
        --output $prefix \\
        --fasta $fasta_name

    EX_COLS=\$(head -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    OBS_COLS=\$(tail -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    if [ "\$EX_COLS" != "\$OBS_COLS" ]; then
        sed -i 's/NA\$/NA\\tNA/' ${prefix}_table.tsv
    fi

    # Cleanup
    rm -rf ${fasta_name} ${fasta_name}.ndb ${fasta_name}.not ${fasta_name}.ntf ${fasta_name}.nto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcroni: $VERSION
    END_VERSIONS
    """
}
