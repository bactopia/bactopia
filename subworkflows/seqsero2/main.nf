/**
 * Predict Salmonella serotypes from genome assemblies.
 *
 * This subworkflow uses [SeqSero2](https://github.com/denglab/SeqSero2) to predict
 * the serotypes of *Salmonella* strains from assembled genomes. It processes each
 * sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords salmonella, serotype, prediction, foodborne, enteric
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation seqsero2
 *
 * @modules csvtk_concat, seqsero2
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: SeqSero2 serotype prediction results in TSV format
 * - `txt`: SeqSero2 serotype prediction results in text format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SEQSERO2 as SEQSERO2_MODULE } from '../../modules/seqsero2/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                 } from 'plugin/nf-bactopia'

workflow SEQSERO2 {
    take:
    seqs: Channel<Record>

    main:
    ch_seqsero2 = SEQSERO2_MODULE(seqs)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_seqsero2, 'tsv', [name: 'seqsero2']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_seqsero2
    run_outputs = ch_csvtk_concat
}
