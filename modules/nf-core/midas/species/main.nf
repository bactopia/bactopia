// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'sra-human-scrubber')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::midas=1.3.2=pyh7cba7a3_7"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process MIDAS_SPECIES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/midas:1.3.2--pyh7cba7a3_7' :
        'quay.io/biocontainers/midas:1.3.2--pyh7cba7a3_7' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}.midas.tsv"), emit: tsv
    tuple val(meta), path("*.abundances.txt")   , emit: abundances
    path "*.{log,err}"                          , emit: logs, optional: true
    path ".command.*"                           , emit: nf_logs
    path "versions.yml"                         , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def read_opts = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        MIDAS_DB=\$(find database/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    else
        MIDAS_DB=\$(find $db/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    fi

    env
    run_midas.py \\
        species \\
        results \\
        $read_opts \\
        $options.args \\
        -d \${MIDAS_DB} \\
        -t $task.cpus

    mv results/species/species_profile.txt ${prefix}.midas.abundances.txt
    midas-summary.py ${prefix} ${prefix}.midas.abundances.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: $VERSION
    END_VERSIONS
    """
}
