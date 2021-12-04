// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'ismapper')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process ISMAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir: query_base) }

    conda (params.enable_conda ? "bioconda::ismapper=2.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ismapper:2.0.2--pyhdfd78af_1' :
        'quay.io/biocontainers/ismapper:2.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    path(reference)
    path(query)

    output:
    tuple val(meta), path("${query_base}/*"), emit: results
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
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
        --output_dir results/ \\
        --queries $query_name \\
        --log ${prefix} \\
        --reference $reference_name \\
        --reads $reads

    mkdir ${query_base}
    mv results/${meta.id}/* ${query_base}

    cat <<-END_VERSIONS > versions.yml
    ismapper:
        ismapper: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
