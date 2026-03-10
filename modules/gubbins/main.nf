/**
 * Detect recombination and construct a recombination-free phylogeny.
 *
 * Uses [Gubbins](https://github.com/nickjcroucher/gubbins) (Genealogies Unbiased By recomBinations In Nucleotide Sequences)
 * to iteratively identify and mask recombinant regions in a multiple sequence alignment. It generates
 * a phylogenetic tree based only on vertically inherited point mutations.
 *
 * @status stable
 * @keywords bacteria, recombination, phylogeny, alignment, msa, evolution, snp
 * @tags complexity:complex input-type:single output-type:multiple features:compression
 * @citation gubbins
 *
 * @input tuple(meta, msa)
 * - `meta`: Groovy Map containing sample information
 * - `msa`: Multiple sequence alignment in FASTA format
 *
 * @output record(meta, masked_aln, fasta, gff, vcf, stats, phylip, embl_predicted, embl_branch, tree, tree_labelled, bootstrap_tree, results, logs, nf_logs, versions)
 * - `masked_aln`: The input alignment with recombinant regions masked (*.masked.aln.gz)
 * - `fasta`: Gubbins internal FASTA alignment
 * - `gff`: Predictions of recombination events in GFF format
 * - `vcf`: SNP calls in VCF format
 * - `stats`: Per-branch statistics on recombination events
 * - `phylip`: Alignment in PHYLIP format
 * - `embl_predicted`: Recombination predictions in EMBL format
 * - `embl_branch`: Branch base reconstruction in EMBL format
 * - `tree`: The final phylogenetic tree constructed from point mutations (*.final_tree.tre)
 * - `tree_labelled`: The final tree with internal node labels
 * - `bootstrap_tree`: The final tree with bootstrap support values (optional)
 */
nextflow.preview.types = true

process GUBBINS {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, msa: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        masked_aln: file("${prefix}.masked.aln.gz"),
        fasta: file("gubbins/${prefix}.filtered_polymorphic_sites.fasta.gz"),
        gff: file("gubbins/${prefix}.recombination_predictions.gff.gz"),
        vcf: file("gubbins/${prefix}.summary_of_snp_distribution.vcf.gz"),
        stats: file("gubbins/${prefix}.per_branch_statistics.csv"),
        phylip: file("gubbins/${prefix}.filtered_polymorphic_sites.phylip"),
        embl_predicted: file("gubbins/${prefix}.recombination_predictions.embl.gz"),
        embl_branch: file("gubbins/${prefix}.branch_base_reconstruction.embl.gz"),
        tree: file("gubbins/${prefix}.final_tree.tre"),
        tree_labelled: file("gubbins/${prefix}.node_labelled.final_tree.tre"),
        bootstrap_tree: file("gubbins/${prefix}.final_bootstrapped_tree.tre", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.masked.aln.gz"),
            files("gubbins/${prefix}.filtered_polymorphic_sites.fasta.gz"),
            files("gubbins/${prefix}.recombination_predictions.gff.gz"),
            files("gubbins/${prefix}.summary_of_snp_distribution.vcf.gz"),
            files("gubbins/${prefix}.per_branch_statistics.csv"),
            files("gubbins/${prefix}.filtered_polymorphic_sites.phylip"),
            files("gubbins/${prefix}.recombination_predictions.embl.gz"),
            files("gubbins/${prefix}.branch_base_reconstruction.embl.gz"),
            files("gubbins/${prefix}.final_tree.tre"),
            files("gubbins/${prefix}.node_labelled.final_tree.tre"),
            files("gubbins/${prefix}.final_bootstrapped_tree.tre", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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
