// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'gtdb')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process GTDBTK_SETUPDB {
    tag "${meta.assembler}-${meta.id}"
    label 'process_high'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::gtdbtk=1.7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:1.7.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gtdbtk:1.7.0--pyhdfd78af_0' }"

    input:
    env GTDBTK_DATA_PATH from params.gtdb

    output:
    path "gtdb-setup.txt"                   , emit: setup
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    """
    if [ "${params.download_gtdb}" == "true" ]; then
        rm -rf !{params.gtdb}/*.tar.gz
        download-db.sh
    else
        echo "skipping GTDB database download"
    fi

    if [ "${params.skip_check}" == "false" ]; then
        gtdbtk check_install && touch gtdb-setup.txt
    else
        echo "skipping GTDB database checks"
        touch gtdb-setup.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    gtdbtk_setupdb:
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
