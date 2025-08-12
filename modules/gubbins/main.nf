process GUBBINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

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
    tuple val(meta), path("*.final_bootstrapped_tree.tre")       , emit: bootstrap_tree, optional: true
    tuple val(meta), path("*.{log,err}")                         , emit: logs, optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = msa.getName().endsWith(".gz") ? true : false
    def msa_name = msa.getName().replace(".gz", "")
    """
    echo "task.ext.args: ${task.ext.args}"
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $msa > $msa_name
    fi

    run_gubbins.py \\
        --threads $task.cpus \\
        --prefix $prefix \\
        $args \\
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
