// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'bakta')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::bakta=1.4.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BAKTA_RUN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.4.0--pyhdfd78af_1' :
        'quay.io/biocontainers/bakta:1.4.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path replicons

    output:
    tuple val(meta), path("${prefix}.{ffn,ffn.gz"), path("${prefix}.{faa,faa.gz"), emit: annotations
    tuple val(meta), path("${prefix}.{embl,embl.gz}")            , emit: embl
    tuple val(meta), path("${prefix}.{faa,faa.gz")               , emit: faa
    tuple val(meta), path("${prefix}.{ffn,ffn.gz")               , emit: ffn
    tuple val(meta), path("${prefix}.{fna,fna.gz}")              , emit: fna
    tuple val(meta), path("${prefix}.{gbff,gbff.gz}")            , emit: gbff
    tuple val(meta), path("${prefix}.{gff3,gff3.gz}")            , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv")         , emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.{faa,faa.gz}"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")                       , emit: tsv
    tuple val(meta), path("${prefix}.txt")                       , emit: txt
    path "*.{log,err}"                                           , emit: logs, optional: true
    path ".command.*"                                            , emit: nf_logs
    path "versions.yml"                                          , emit: versions


    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    def replicons_opt = replicons ? "--replicons ${replicons[0]}" : ""
    """
    bakta \\
        $options.args \\
        --threads $task.cpus \\
        --prefix ${prefix} \\
        --db $db \\
        $proteins_opt \\
        $prodigal_tf \\
        $replicons_opt \\
        $fasta

    if [[ "!{params.skip_compression}" == "false" ]]; then
        gzip --best ${prefix}.embl
        gzip --best ${prefix}.faa
        gzip --best ${prefix}.ffn
        gzip --best ${prefix}.fna
        gzip --best ${prefix}.gbff
        gzip --best ${prefix}.gff3
        gzip --best ${prefix}.hypotheticals.faa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
