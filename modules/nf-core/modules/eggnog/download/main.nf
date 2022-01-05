// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'eggnog_download')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process EGGNOG_DOWNLOAD {
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.6--pyhdfd78af_0' :
        'quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0' }"

    output:
    path("eggnog/*")                        , emit: db
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs    , optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    """
    mkdir eggnog
    download_eggnog_data.py \\
        $options.args \\
        -y \\
        --data_dir eggnog/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """
}
