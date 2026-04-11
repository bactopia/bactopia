/**
 * Build a pangenome from GFF3 annotations using PIRATE.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [PIRATE](https://github.com/SionBayliss/PIRATE).
 * PIRATE is a scalable pangenome toolbox that clusters orthologous genes at multiple identity thresholds.
 * It is particularly useful for highly diverse datasets as it can handle divergent gene families
 * and provides flexible clustering options for different analytical needs.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation pirate
 *
 * @modules pirate as pirate_module
 *
 * @input record(meta, gff)
 * - `meta`: Groovy Record containing sample information
 * - `gff`: Set of GFF3 annotation files representing the genomic annotations for each sample
 *
 * @output sample_outputs
 *   - `aln`: Core genome alignment in FASTA format (optional)
 *   - `csv`: Gene presence/absence matrix in CSV format
 *   - `supplemental`: Directory containing PIRATE intermediate files and detailed outputs
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'

workflow PIRATE {
    take:
    gff: Channel<Record>

    main:
    ch_pirate = PIRATE_MODULE(gff)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_pirate
}
