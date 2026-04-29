/**
 * Run all available MLST schemes for a species against an assembly
 *
 * This subworkflow uses [GigaTyper](https://github.com/rpetit3/gigatyper) to run all available mlst schemes for a species against an assembly.
 * It processes each sample individually and aggregates the results into
 * a single consolidated report.
 *
 * @status stable
 * @keywords mlst, typing, multi-scheme
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,
 * @citation gigatyper
 *
 * @modules csvtk_concat, gigatyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: MLST results across all schemes
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with gigatyper results from all samples
 */
nextflow.enable.types = true

include { GIGATYPER as GIGATYPER_MODULE } from '../../modules/gigatyper/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow GIGATYPER {
    take:
    fna: Channel<Record>

    main:
    ch_gigatyper = GIGATYPER_MODULE(fna)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_gigatyper, 'tsv', [name: 'gigatyper']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_gigatyper
    run_outputs = ch_csvtk_concat
}
