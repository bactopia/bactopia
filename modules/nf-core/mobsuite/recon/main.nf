// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'mobsuite')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::mob_suite=3.1.8"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process MOBSUITE_RECON {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.8--pyhdfd78af_1' :
        'quay.io/biocontainers/mob_suite:3.1.8--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/chromosome.fasta") , emit: chromosome
    tuple val(meta), path("results/contig_report.txt"), emit: contig_report
    tuple val(meta), path("results/plasmid_*.fasta")  , emit: plasmids        , optional: true
    tuple val(meta), path("results/${prefix}-mobtyper.txt") , emit: mobtyper_results, optional: true
    path "*.{log,err}"                                , emit: logs, optional: true
    path ".command.*"                                 , emit: nf_logs
    path "versions.yml"                               , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        $options.args \\
        --num_threads $task.cpus \\
        --outdir results \\
        --sample_id $prefix

    if [[ -f "results/mobtyper_results.txt" ]]; then
        mv results/mobtyper_results.txt results/${prefix}-mobtyper.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
