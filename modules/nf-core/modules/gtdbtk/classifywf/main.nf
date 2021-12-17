// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'gtdb')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::gtdbtk=1.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:1.7.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gtdbtk:1.7.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fna, stageAs: 'fna-tmp/*')
    path db, stageAs: 'gtdb/*'

    output:
    path "results/*"                        , emit: results
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    export GTDBTK_DATA_PATH="${params.gtdb}"
    mkdir fna
    cp -L fna-tmp/* fna/
    find fna/ -name "*.fna.gz" | xargs gunzip

    gtdbtk classify_wf \\
        $options.args \\
        --cpus $task.cpus \\
        --pplacer_cpus $task.cpus \\
        --genome_dir ./fna \\
        --out_dir results \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    gtdbtk_classifywf:
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
