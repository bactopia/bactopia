// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'pneumocat')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::pneumocat=1.2.1 bioconda::bowtie2==2.3.5"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

process PNEUMOCAT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pneumocat:1.2.1--0' :
        'quay.io/biocontainers/pneumocat:1.2.1--0' }"

    input:
    tuple val(meta), path(reads)

    when:
    meta.single_end == false

    output:
    tuple val(meta), path("*.xml")               , optional: true, emit: xml
    tuple val(meta), path("coverage_summary.txt"), optional: true, emit: txt
    path "*.{stdout,stderr}"                     , optional: true, emit: logs
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        --threads $task.cpus \\
        --output_dir ./

    # clean up

    # PneumoCAT uses first match in a glob, so moves between R1 and R2
    if [ -f ${prefix}_R1.results.xml ]; then
        mv ${prefix}_R1.results.xml ${prefix}.results.xml
    else
        mv ${prefix}_R2.results.xml ${prefix}.results.xml
    fi
    mv logs/* ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """
}
