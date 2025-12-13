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
 * @output masked_aln      The input alignment with recombinant regions masked (*.masked.aln.gz)
 * @output fasta           Gubbins internal FASTA alignment
 * @output gff             Predictions of recombination events in GFF format
 * @output vcf             SNP calls in VCF format
 * @output stats           Per-branch statistics on recombination events
 * @output phylip          Alignment in PHYLIP format
 * @output embl_predicted  Recombination predictions in EMBL format
 * @output embl_branch     Branch base reconstruction in EMBL format
 * @output tree            The final phylogenetic tree constructed from point mutations (*.final_tree.tre)
 * @output tree_labelled   The final tree with internal node labels
 * @output bootstrap_tree  The final tree with bootstrap support values (optional)
 * @output logs            Optional software execution logs containing warnings/errors
 * @output nf_logs         Nextflow execution scripts and logs for debugging
 * @output versions        A YAML formatted file with software versions
 */
nextflow.preview.types = true

process GUBBINS {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, msa) : Tuple<Map, Set<Path>>

    output:
    masked_aln     = tuple(meta, files("*.masked.aln.gz"))
    fasta          = tuple(meta, files("gubbins/*.fasta.gz"))
    gff            = tuple(meta, files("gubbins/*.gff.gz"))
    vcf            = tuple(meta, files("gubbins/*.vcf.gz"))
    stats          = tuple(meta, files("gubbins/*.csv"))
    phylip         = tuple(meta, files("gubbins/*.phylip"))
    embl_predicted = tuple(meta, files("gubbins/*.recombination_predictions.embl.gz"))
    embl_branch    = tuple(meta, files("gubbins/*.branch_base_reconstruction.embl.gz"))
    tree           = tuple(meta, files("gubbins/*.final_tree.tre"))
    tree_labelled  = tuple(meta, files("gubbins/*.node_labelled.final_tree.tre"))
    bootstrap_tree = tuple(meta, files("gubbins/*.final_bootstrapped_tree.tre", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, files("versions.yml"))

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
    def is_compressed = msa.toList()[0].getName().endsWith(".gz") ? true : false
    def msa_name = msa.toList()[0].getName().replace(".gz", "")
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
