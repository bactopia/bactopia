// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'hpsuissero')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::hpsuissero=1.0.1" 
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
def VERSION   = '1.0.1' // Version information not provided by tool on CLI

process HPSUISSERO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hpsuissero:1.0.1--hdfd78af_0' :
        'quay.io/biocontainers/hpsuissero:1.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    HpsuisSero.sh \\
        -i $fasta_name \\
        -o ./ \\
        -s $prefix \\
        -x fasta \\
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpsuissero: $VERSION
    END_VERSIONS
    """
}
