// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'gubbins')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::gubbins=3.3.4"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process GUBBINS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.3.4--py310pl5321he4a0461_0' :
        'quay.io/biocontainers/gubbins:3.3.4--py310pl5321he4a0461_0' }"

    input:
    tuple val(meta), path(msa)

    output:
    tuple val(meta), path("*.masked.aln.gz")                     , emit: masked_aln
    tuple val(meta), path("*.fasta.gz")                          , emit: fasta
    tuple val(meta), path("*.gff.gz")                            , emit: gff
    tuple val(meta), path("*.vcf.gz")                            , emit: vcf
    tuple val(meta), path("*.csv")                               , emit: stats
    tuple val(meta), path("*.phylip")                            , emit: phylip
    tuple val(meta), path("*.recombination_predictions.embl.gz") , emit: embl_predicted
    tuple val(meta), path("*.branch_base_reconstruction.embl.gz"), emit: embl_branch
    tuple val(meta), path("*.final_tree.tre")                    , emit: tree
    tuple val(meta), path("*.node_labelled.final_tree.tre")      , emit: tree_labelled
    path "*.{log,err}"                           , optional: true, emit: logs
    path ".command.*"                                            , emit: nf_logs
    path "versions.yml"                                          , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    run_gubbins.py \\
        --threads $task.cpus \\
        --prefix $prefix \\
        $options.args \\
        $msa_name

    # Create masked alignment
    mask_gubbins_aln.py \\
        --aln $msa_name \\
        --gff ${prefix}.recombination_predictions.gff \\
        --out ${prefix}.masked.aln

    # Cleanup
    gzip *.masked.aln *.embl *.fasta *.gff *.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}
