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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation blast
 *
 * @modules tblastn, csvtk_concat
 *
 * @input record(meta, blastdb)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
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
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTN {
    take:
    assembly: Channel<Record>
    query: Path

    main:
    TBLASTN_MODULE(assembly, query)
    CSVTK_CONCAT(gather(TBLASTN_MODULE.out, 'tblastn', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = TBLASTN_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
