/**
 * Identify and type Stx operons from assembled genomic sequences
 *
 * This subworkflow uses [StxTyper](https://github.com/ncbi/stxtyper) to identify and type stx operons from assembled genomic sequences.
 * It processes each sample individually and aggregates the results into
 * a single consolidated report.
 *
 * @status stable
 * @keywords stx, shiga toxin, typing, stec, virulence
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation stxtyper
 *
 * @modules csvtk_concat, stxtyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited Stx operon typing results
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with stxtyper results from all samples
 */
nextflow.enable.types = true

include { STXTYPER as STXTYPER_MODULE } from '../../modules/stxtyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                 } from 'plugin/nf-bactopia'

workflow STXTYPER {
    take:
    fna: Channel<Record>

    main:
    ch_stxtyper = STXTYPER_MODULE(fna)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_stxtyper, 'tsv', [name: 'stxtyper']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_stxtyper
    run_outputs = ch_csvtk_concat
}
