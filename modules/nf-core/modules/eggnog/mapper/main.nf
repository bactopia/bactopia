// Import generic module functions
include { initOptions; saveFiles } from '../../../../../lib/nf/functions'

options     = initOptions(params.options ? params.options : [:], 'eggnog_mapper')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process EGGNOG_MAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.6--pyhdfd78af_0' :
        'quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0' }"

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
        --data_dir ./ \\
        --output $prefix \\
        -i $fasta

    cat <<-END_VERSIONS > versions.yml
    eggnog_mapper:
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """
}
