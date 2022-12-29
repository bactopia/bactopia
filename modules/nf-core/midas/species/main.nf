// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'midas')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::midas=1.3.2"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process MIDAS_SPECIES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/midas:1.3.2--pyh7cba7a3_7' :
        'quay.io/biocontainers/midas:1.3.2--pyh7cba7a3_7' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
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

    run_midas.py \\
        species \\
        results \\
        $read_opts \\
        $options.args \\
        -d \$MIDAS_DB \\
        -t $task.cpus

    mv results/species/species_profile.txt ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: $VERSION
    END_VERSIONS
    """
}
