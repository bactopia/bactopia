// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'ismapper')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::ismapper=2.0.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ISMAPPER {
    tag "$prefix"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ismapper:2.0.2--pyhdfd78af_1' :
        'quay.io/biocontainers/ismapper:2.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    path(reference)
    path(query)

    output:
    tuple val(meta), path("results/*"), emit: results
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def ref_compressed = reference.getName().endsWith(".gz") ? true : false
    def reference_name = reference.getName().replace(".gz", "")
    def query_compressed = query.getName().endsWith(".gz") ? true : false
    def query_name = query.getName().replace(".gz", "")
    query_base = query.getSimpleName()
    """
    if [ "$ref_compressed" == "true" ]; then
        gzip -c -d $reference > $reference_name
    fi
    if [ "$query_compressed" == "true" ]; then
        gzip -c -d $query > $query_name
    fi
    
    ismap \\
        $options.args \\
        --t $task.cpus \\
        --output_dir $prefix \\
        --queries $query_name \\
        --log ${prefix} \\
        --reference $reference_name \\
        --reads $reads

    # Reorganize output files
    mkdir results
    mv $prefix/*/* results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ismapper: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
