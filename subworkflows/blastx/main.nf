/**
 * Translate nucleotide sequences and search protein database.
 *
 * This subworkflow uses [BLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastx)
 * from the NCBI BLAST+ suite to translate nucleotide sequences in all six reading frames
 * and search them against a protein database. It processes each assembly individually
 * and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords blast, protein, translation, alignment, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation blast
 *
 * @modules blastx, csvtk_concat
 *
 * @input record(meta, blastdb)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `blastdb`: A compressed tarball containing the nucleotide BLAST database
 *
 * @input query
 * Path to BLAST protein database for searching translated sequences
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of alignments (standard BLAST outfmt 6)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                        } from 'plugin/nf-bactopia'

workflow BLASTX {
    take:
    blastdb: Channel<Record>
    query: Path

    main:
    BLASTX_MODULE(blastdb, query)
    CSVTK_CONCAT(gatherCsvtk(BLASTX_MODULE.out, 'tsv', [name: 'blastx']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = BLASTX_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
