/**
 * Rapidly identify genes by creating local assemblies from paired-end reads.
 *
 * This subworkflow uses [ARIBA](https://github.com/sanger-pathogens/ariba)
 * (Antimicrobial Resistance Identification By Assembly) to rapidly identify genes
 * in a database by creating local assemblies. It first downloads and prepares an ARIBA database,
 * then analyzes paired-end reads to identify genes, and finally aggregates results across all samples.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords bacteria, reads, antimicrobial resistance, virulence, local assembly
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent, resource-download
 * @citation ariba
 *
 * @modules ariba_getref, ariba_run, csvtk_concat
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by ARIBA)
 * - `lr`: Long reads (not supported by ARIBA)
 *
 * @input db
 * Database name for ARIBA analysis (e.g., ncbi, card, vfdb, resfinder, argannot)
 *
 * @output sample_outputs
 * - `report`: Tab-delimited detailed report of gene detection results
 * - `summary`: Comma-separated condensed summary of detected genes
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { ARIBA_GETREF                         } from '../../modules/ariba/getref/main'
include { ARIBA_RUN                            } from '../../modules/ariba/run/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_REPORT  } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SUMMARY } from '../../modules/csvtk/concat/main'
include { gather                               } from 'plugin/nf-bactopia'

workflow ARIBA {
    take:
    reads: Channel<Record>
    db: String

    main:
    ARIBA_GETREF(db)
    ARIBA_RUN(reads, ARIBA_GETREF.out.map { r -> r.db })
    CSVTK_CONCAT_REPORT(gather(ARIBA_RUN.out, 'report', [name: "${db}-report", args: '-C "$" --lazy-quotes']), 'tsv', 'tsv')
    CSVTK_CONCAT_SUMMARY(gather(ARIBA_RUN.out, 'summary', [name: "${db}-summary", args: '--lazy-quotes']), 'csv', 'csv')

    emit:
    // Published outputs
    sample_outputs = ARIBA_RUN.out
    run_outputs = CSVTK_CONCAT_REPORT.out.mix(CSVTK_CONCAT_SUMMARY.out)
}
