/**
 * Search protein query sequences against nucleotide database.
 *
 * This subworkflow uses [TBLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=tblastn)
 * from the NCBI BLAST+ suite to search protein query sequences against a nucleotide database
 * translated in all six reading frames. It processes each assembly individually
 * and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords blast, protein, nucleotide, alignment, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation blast
 *
 * @modules blast_tblastn, csvtk_concat
 *
 * @input record(meta, blastdb)
 * - `meta`: Groovy Record containing sample information
 * - `blastdb`: A compressed tarball containing the nucleotide BLAST database
 *
 * @input query
 * Path to protein query sequences for searching against translated nucleotide database
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited protein-to-translated nucleotide alignment results (BLAST outfmt 6)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../modules/blast/tblastn/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                     } from 'plugin/nf-bactopia'

workflow TBLASTN {
    take:
    blastdb: Channel<Record>
    query: Path

    main:
    ch_tblastn = TBLASTN_MODULE(blastdb, query)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_tblastn, 'tsv', [name: 'tblastn']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_tblastn
    run_outputs = ch_csvtk_concat
}
