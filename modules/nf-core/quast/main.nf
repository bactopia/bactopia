// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'quast')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::quast=5.2.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_2' :
        'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_2' }"

    input:
    tuple val(meta), path(fasta), path(meta_file)

    output:
    tuple val(meta), path("results/${meta.id}.tsv"), emit: tsv
    tuple val(meta), path("results/*")             , emit: results
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    est_ref_size=\$(tail -n 1 $meta_file | cut -f 7)
    if [ "\${est_ref_size}" != "0" ]; then
        est_ref_size="--est-ref-size \${est_ref_size}"
    fi

    quast ${fasta_name} \${est_ref_size} \\
        -o results \\
        --threads ${task.cpus} \\
        $options.args \\
        --glimmer

    mv results/quast.log ./
    mv results/transposed_report.tsv results/${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
