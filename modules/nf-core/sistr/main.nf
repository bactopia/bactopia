// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'sistr')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::sistr_cmd=1.1.1"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SISTR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sistr_cmd:1.1.1--pyh864c0ab_2' :
        'quay.io/biocontainers/sistr_cmd:1.1.1--pyh864c0ab_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")            , emit: tsv
    tuple val(meta), path("*-allele.fasta.gz"), emit: allele_fasta
    tuple val(meta), path("*-allele.json.gz") , emit: allele_json
    tuple val(meta), path("*-cgmlst.csv")     , emit: cgmlst_csv
    path "*.{log,err}"                        , emit: logs, optional: true
    path ".command.*"                         , emit: nf_logs
    path "versions.yml"                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    sistr \\
        --qc \\
        $options.args \\
        --threads $task.cpus \\
        --alleles-output ${prefix}-allele.json \\
        --novel-alleles ${prefix}-allele.fasta \\
        --cgmlst-profiles ${prefix}-cgmlst.csv \\
        --output-prediction ${prefix} \\
        --output-format tab \\
        $fasta_name

    mv ${prefix}.tab ${prefix}.tsv
    gzip ${prefix}-allele.json
    gzip ${prefix}-allele.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sistr: \$(echo \$(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' )
    END_VERSIONS
    """
}
