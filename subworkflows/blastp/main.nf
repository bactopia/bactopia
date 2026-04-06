/**
 * Search protein sequences against protein database.
 *
 * This subworkflow uses [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastp)
 * from the NCBI BLAST+ suite to search protein sequences against a protein database.
 * It processes each assembly individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords blast, protein, alignment, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation blast
 *
 * @modules blast_blastp, csvtk_concat
 *
 * @input record(meta, blastdb)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `blastdb`: A compressed tarball containing the protein BLAST database
 *
 * @input query
 * Path to protein database for searching against translated sequences
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of alignments (standard BLAST outfmt 6)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BLAST_BLASTP as BLASTP_MODULE } from '../../modules/blast/blastp/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow BLASTP {
    take:
    blastdb: Channel<Record>
    query: Value<Path>

    main:
    BLASTP_MODULE(blastdb, query)
    CSVTK_CONCAT(gatherCsvtk(BLASTP_MODULE.out, 'tsv', [name: 'blastp']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = BLASTP_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
