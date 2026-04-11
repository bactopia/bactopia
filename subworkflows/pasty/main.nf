/**
 * Predict serogroups of Pseudomonas aeruginosa from assemblies.
 *
 * This subworkflow uses [Pasty](https://github.com/rpetit3/pasty) to perform in silico
 * serogrouping of *Pseudomonas aeruginosa* isolates from assembled genomes. It identifies
 * O-antigen biosynthesis genes to classify isolates into their known serogroups using
 * BLAST-based homology searches.
 *
 * @status stable
 * @keywords pseudomonas aeruginosa, serogroup, typing, o-antigen, prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation pasty
 *
 * @modules pasty, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary file with the predicted O-antigen serogroup
 * - `blast`: A tab-delimited file of all raw BLAST hits used for the prediction
 * - `details`: A tab-delimited file with detailed gene hits for each serogroup tested
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk           } from 'plugin/nf-bactopia'

workflow PASTY {
    take:
    assembly: Channel<Record>

    main:
    ch_pasty = PASTY_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_pasty, 'tsv', [name: 'pasty']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_pasty
    run_outputs = ch_csvtk_concat
}
