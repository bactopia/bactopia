// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process EGGNOG_MAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}",
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
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.emapper.hits")                , emit: hits
    tuple val(meta), path("*.emapper.seed_orthologs")      , emit: seed_orthologs
    tuple val(meta), path("*.emapper.annotations")         , emit: annotations
    tuple val(meta), path("*.emapper.annotations.xlsx")    , emit: xlsx     , optional: true
    tuple val(meta), path("*.emapper.orthologs")           , emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.genepred.fasta")      , emit: genepred , optional: true
    tuple val(meta), path("*.emapper.gff")                 , emit: gff      , optional: true
    tuple val(meta), path("*.emapper.no_annotations.fasta"), emit: no_anno  , optional: true
    tuple val(meta), path("*.emapper.pfam")                , emit: pfam     , optional: true
    path "*.{stdout.txt,stderr.txt,log,err}"               , emit: logs     , optional: true
    path ".command.*"                                      , emit: nf_logs
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    emapper.py \\
        $options.args \\
        --cpu $task.cpus \\
        --data_dir $db \\
        --output $prefix \\
        -i $fasta

    cat <<-END_VERSIONS > versions.yml
    eggnog_mapper:
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """

    stub:
    """
    touch test.emapper.hits
    touch test.emapper.seed_orthologs
    touch test.emapper.annotations

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """
}
