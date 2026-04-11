/**
 * Build a pangenome from GFF3 annotations using Panaroo.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [Panaroo](https://github.com/gtonkinhill/panaroo).
 * Panaroo is a pangenome pipeline that produces polished pangenomes by removing errors and
 * contamination from input annotations. It generates gene presence/absence matrices and core-genome
 * alignments suitable for downstream phylogenetic analysis.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation panaroo
 *
 * @modules panaroo_run
 *
 * @input record(meta, gff)
 * - `meta`: Groovy Record containing sample information
 * - `gff`: Set of GFF3 annotation files representing the genomic annotations for each sample
 *
 * @output sample_outputs
 *   - `aln`: Core genome alignment in FASTA format (optional)
 *   - `filtered_aln`: Core genome alignment with recombinant regions filtered out (optional)
 *   - `csv`: Gene presence/absence matrix in Roary-compatible CSV format (optional)
 *   - `panaroo_csv`: Gene presence/absence matrix in Panaroo's native CSV format (optional)
 *   - `supplemental`: Directory containing Panaroo intermediate files and data structures
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { PANAROO_RUN } from '../../modules/panaroo/run/main'

workflow PANAROO {
    take:
    gff: Channel<Record>

    main:
    ch_panaroo_run = PANAROO_RUN(gff)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_panaroo_run
}
