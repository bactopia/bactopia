// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'plasmidfinder')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::plasmidfinder=2.1.6=py310hdfd78af_1"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
def VERSION = '2.1.6' // Version information not provided by tool on CLI

process PLASMIDFINDER {
    tag "$meta.id"
    label 'process_low'

    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidfinder:2.1.6--py310hdfd78af_1' :
        'quay.io/biocontainers/plasmidfinder:2.1.6--py310hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json")                 , emit: json
    tuple val(meta), path("*.txt")                  , emit: txt
    tuple val(meta), path("*.tsv")                  , emit: tsv
    tuple val(meta), path("*-hit_in_genome_seq.fsa"), emit: genome_seq
    tuple val(meta), path("*-plasmid_seqs.fsa")     , emit: plasmid_seq
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    plasmidfinder.py \\
        $options.args \\
        -i $fasta_name \\
        -o ./ \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv data.json ${prefix}.json
    mv results.txt ${prefix}.txt
    mv results_tab.tsv ${prefix}.tsv
    mv Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${prefix}-plasmid_seqs.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: $VERSION
    END_VERSIONS
    """
}
