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
 * @input record(meta, aln)
 * - `meta`: Groovy Map containing sample information
 * - `aln`: Multiple sequence alignment in FASTA format
 *
 * @output record(meta, masked_aln, results, logs, nf_logs, versions)
 * - `masked_aln`: The input alignment with recombinant regions masked (*.masked.aln.gz)
 *
 * @results gubbins
 * - `*.branch_base_reconstruction.embl.gz`: Per-branch base reconstruction in EMBL format
 * - `*.recombination_predictions.embl.gz`: Predicted recombination events in EMBL format
 * - `*.recombination_predictions.gff.gz`: Predicted recombination regions in GFF format
 * - `*.per_branch_statistics.csv`: Statistics for each branch of the phylogeny
 * - `*.filtered_polymorphic_sites.fasta.gz`: Alignment of filtered polymorphic sites
 * - `*.filtered_polymorphic_sites.phylip`: Filtered polymorphic sites in PHYLIP format
 * - `*.final_tree.tre`: Final recombination-free phylogenetic tree
 * - `*.node_labelled.final_tree.tre`: Final tree with internal node labels
 * - `*.summary_of_snp_distribution.vcf.gz`: SNP distribution summary in VCF format
 */
nextflow.preview.types = true

process GUBBINS {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        aln: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        masked_aln: file("${prefix}.masked.aln.gz"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.masked.aln.gz"),
            files("gubbins/*"),
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = ""
    meta.logs_dir = "${meta.process_name}/logs/"

    def is_compressed = aln.getName().endsWith(".gz") ? true : false
    def aln_name = aln.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${aln} > ${aln_name}
    fi

    run_gubbins.py \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${task.ext.args} \\
        ${aln_name}

    # Create masked alignment
    mask_gubbins_aln.py \\
        --aln ${aln_name} \\
        --gff ${prefix}.recombination_predictions.gff \\
        --out ${prefix}.masked.aln

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${aln_name}
    fi
    gzip *.masked.aln *.embl *.fasta *.gff *.vcf

    # Move supplemental outputs to gubbins folder
    mkdir gubbins
    mv ${prefix}.branch_base_reconstruction.embl.gz \
       ${prefix}.recombination_predictions.embl.gz \
       ${prefix}.recombination_predictions.gff.gz \
       ${prefix}.per_branch_statistics.csv \
       ${prefix}.filtered_polymorphic_sites.fasta.gz \
       ${prefix}.filtered_polymorphic_sites.phylip \
       ${prefix}.final_tree.tre \
       ${prefix}.node_labelled.final_tree.tre \
       ${prefix}.summary_of_snp_distribution.vcf.gz \
       gubbins/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}
