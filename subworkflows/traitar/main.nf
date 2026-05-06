/**
 * Predict phenotypic traits from microbial genomes
 *
 * This subworkflow uses [Traitar](https://github.com/nick-youngblut/traitar3/) to predict phenotypic traits from microbial genomes.
 * It can download and prepare the Pfam database on-demand or use a pre-existing database.
 * It processes each sample individually and aggregates the results into
 * a single consolidated report.
 *
 * @status stable
 * @keywords phenotype, traits, pfam
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent,conditional-logic,resource-download
 * @citation traitar
 *
 * @modules csvtk_concat, traitar_download, traitar_run
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input database
 * Optional pre-existing Pfam-A HMM file
 *
 * @input download_traitar
 * Boolean flag to trigger automatic database download
 *
 * @output sample_outputs
 * - `majority_tsv`: Majority-vote combined phenotype trait predictions
 * - `single_tsv`: Single-votes combined phenotype trait predictions
 *
 * @output run_outputs
 * - `csv`: Merged TSV files with traitar majority-vote and single-vote results from all samples
 */
nextflow.enable.types = true

include { TRAITAR_DOWNLOAD                      } from '../../modules/traitar/download/main'
include { TRAITAR_RUN                           } from '../../modules/traitar/run/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_MAJORITY } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SINGLE   } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                           } from 'plugin/nf-bactopia'

workflow TRAITAR {
    take:
    fna: Channel<Record>
    database: Path?
    download_traitar: Boolean

    main:
    ch_traitar_run = channel.empty()
    if (download_traitar) {
        ch_traitar_download = TRAITAR_DOWNLOAD()
        ch_traitar_run = TRAITAR_RUN(fna, ch_traitar_download.map { r -> r.db })
    } else {
        ch_traitar_run = TRAITAR_RUN(fna, database)
    }
    ch_majority_concat = CSVTK_CONCAT_MAJORITY(gatherCsvtk(ch_traitar_run, 'majority_tsv', [name: 'traitar-majority']), 'tsv', 'tsv')
    ch_single_concat = CSVTK_CONCAT_SINGLE(gatherCsvtk(ch_traitar_run, 'single_tsv', [name: 'traitar-single']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_traitar_run
    run_outputs = ch_majority_concat.mix(ch_single_concat)
}
