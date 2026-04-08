/**
 * Build a pangenome from GFF3 annotations using Roary.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [Roary](https://github.com/sanger-pathogens/Roary).
 * Roary is a rapid pangenome pipeline that processes large numbers of annotated genomes to produce
 * gene presence/absence matrices and core-genome alignments. It is particularly optimized for
 * bacterial datasets and can handle hundreds of genomes efficiently.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation roary
 *
 * @modules roary as roary_module
 *
 * @input record(meta, gff)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `gff`: Set of GFF3 annotation files representing the genomic annotations for each sample
 *
 * @output sample_outputs
 *
 * @output run_outputs
 *   - `aln`: Core genome alignment in FASTA format (optional)
 *   - `csv`: Gene presence/absence table
 *   - `supplemental`: Supplemental files including accessory binary genes and graphs
 */
nextflow.preview.types = true

include { ROARY as ROARY_MODULE } from '../../modules/roary/main'

workflow ROARY {
    take:
    gff: Channel<Record>

    main:
    ch_roary = ROARY_MODULE(gff)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_roary
}
