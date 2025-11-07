process GUBBINS {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(msa)

    output:
    tuple val(meta), path("*.masked.aln.gz")                             , emit: masked_aln
    tuple val(meta), path("gubbins/*.fasta.gz")                          , emit: fasta
    tuple val(meta), path("gubbins/*.gff.gz")                            , emit: gff
    tuple val(meta), path("gubbins/*.vcf.gz")                            , emit: vcf
    tuple val(meta), path("gubbins/*.csv")                               , emit: stats
    tuple val(meta), path("gubbins/*.phylip")                            , emit: phylip
    tuple val(meta), path("gubbins/*.recombination_predictions.embl.gz") , emit: embl_predicted
    tuple val(meta), path("gubbins/*.branch_base_reconstruction.embl.gz"), emit: embl_branch
    tuple val(meta), path("gubbins/*.final_tree.tre")                    , emit: tree
    tuple val(meta), path("gubbins/*.node_labelled.final_tree.tre")      , emit: tree_labelled
    tuple val(meta), path("gubbins/*.final_bootstrapped_tree.tre")       , emit: bootstrap_tree, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = ""
    meta.logs_dir = "${meta.process_name}/logs/"
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    run_gubbins.py \\
        --threads $task.cpus \\
        --prefix $prefix \\
        ${task.ext.args} \\
        $msa_name

    # Create masked alignment
    mask_gubbins_aln.py \\
        --aln $msa_name \\
        --gff ${prefix}.recombination_predictions.gff \\
        --out ${prefix}.masked.aln

    # Cleanup
    gzip *.masked.aln *.embl *.fasta *.gff *.vcf

    # Move outputs to tool specific folder
    mkdir gubbins
    mv ${prefix}* gubbins/
    mv gubbins/${prefix}.masked.aln.gz ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}
