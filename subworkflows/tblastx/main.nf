/**
 * Translate nucleotide query sequences and search nucleotide database.
 *
 * This subworkflow uses [TBLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=tblastx)
 * from the NCBI BLAST+ suite to translate nucleotide query sequences in all six reading frames
 * and search them against a nucleotide database also translated in all six reading frames.
 * It processes each assembly individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords blast, nucleotide, translation, alignment, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation blast
 *
 * @modules blast_tblastx, csvtk_concat
 *
 * @input record(meta, blastdb)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `blastdb`: A compressed tarball containing the nucleotide BLAST database
 *
 * @input query
 * Path to nucleotide query sequences that will be translated and searched against translated database
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited translated nucleotide-to-translated nucleotide alignment results (BLAST outfmt 6)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BLAST_TBLASTX as TBLASTX_MODULE } from '../../modules/blast/tblastx/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                     } from 'plugin/nf-bactopia'

workflow TBLASTX {
    take:
    blastdb: Channel<Record>
    query: Value<Path>

    main:
    ch_tblastx = TBLASTX_MODULE(blastdb, query)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_tblastx, 'tsv', [name: 'tblastx']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_tblastx
    run_outputs = ch_csvtk_concat
}
