nextflow.preview.types = true

process GUBBINS {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, msa) : Tuple<Map, Path>

    output:
    masked_aln     = tuple(meta, file("*.masked.aln.gz"))
    fasta          = tuple(meta, file("gubbins/*.fasta.gz"))
    gff            = tuple(meta, file("gubbins/*.gff.gz"))
    vcf            = tuple(meta, file("gubbins/*.vcf.gz"))
    stats          = tuple(meta, file("gubbins/*.csv"))
    phylip         = tuple(meta, file("gubbins/*.phylip"))
    embl_predicted = tuple(meta, file("gubbins/*.recombination_predictions.embl.gz"))
    embl_branch    = tuple(meta, file("gubbins/*.branch_base_reconstruction.embl.gz"))
    tree           = tuple(meta, file("gubbins/*.final_tree.tre"))
    tree_labelled  = tuple(meta, file("gubbins/*.node_labelled.final_tree.tre"))
    bootstrap_tree = tuple(meta, file("gubbins/*.final_bootstrapped_tree.tre", optional: true))
    logs           = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin       = tuple(meta, file(".command.begin"))
    nf_err         = tuple(meta, file(".command.err"))
    nf_log         = tuple(meta, file(".command.log"))
    nf_out         = tuple(meta, file(".command.out"))
    nf_run         = tuple(meta, file(".command.run"))
    nf_sh          = tuple(meta, file(".command.sh"))
    nf_trace       = tuple(meta, file(".command.trace"))
    versions       = tuple(meta, file("versions.yml"))

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
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${msa} > ${msa_name}
    fi

    run_gubbins.py \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${task.ext.args} \\
        ${msa_name}

    # Create masked alignment
    mask_gubbins_aln.py \\
        --aln ${msa_name} \\
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
