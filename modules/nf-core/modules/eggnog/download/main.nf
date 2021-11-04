// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process EGGNOG_DOWNLOAD {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), subworkflow: options.subworkflow, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.6--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0"
    }

    input:
    val(meta)

    output:
    path("*.db*")                           , emit: db
    path("*.dmnd")                          , emit: proteins, optional: true
    path("hmmer/")                          , emit: hmmer   , optional: true
    path("mmseqs/")                         , emit: mmseqs  , optional: true
    path("pfam/")                           , emit: pfam    , optional: true
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs    , optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    download_eggnog_data.py \\
        $options.args \\
        -y \\
        --data_dir ./

    cat <<-END_VERSIONS > versions.yml
    eggnog_download:
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """

    stub:
    """
    touch eggnog.db
    touch eggnog.taxa.db
    touch eggnog.taxa.db.traverse.pkl

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """
}
