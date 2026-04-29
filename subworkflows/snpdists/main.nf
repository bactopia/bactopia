/**
 * Calculate pairwise SNP distances from sequence alignments.
 *
 * This subworkflow uses [snp-dists](https://github.com/tseemann/snp-dists) to compute
 * pairwise SNP distance matrices from multiple sequence alignments. It reads an
 * alignment file (typically a core-genome alignment) and calculates the number
 * of SNP differences between each pair of sequences, producing a distance matrix
 * useful for phylogenetic and epidemiological analyses.
 *
 * @status stable
 * @keywords snp, distance, alignment, phylogeny, core-genome
 * @tags complexity:simple input-type:single output-type:multiple
 * @citation snpdists
 *
 * @modules snpdists
 *
 * @input record(meta, aln)
 * - `meta`: Groovy Record containing sample information
 * - `aln`: Multiple sequence alignment in FASTA format
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `tsv`: Pairwise SNP distance matrix in TSV format
 */
nextflow.enable.types = true

include { SNPDISTS as SNPDISTS_MODULE } from '../../modules/snpdists/main'

workflow SNPDISTS {
    take:
    alignment: Channel<Record>

    main:
    ch_snpdists = SNPDISTS_MODULE(alignment)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_snpdists
}
