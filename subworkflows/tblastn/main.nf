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
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format that will be translated and searched
 *
 * @input query
 * Path to protein query sequences for searching against translated nucleotide database
 *
 * @output sample_outputs  Record outputs from TBLASTN module per sample
 * @output run_outputs   Record with merged CSV and associated logs/versions
 */
nextflow.preview.types = true

include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../modules/blast/tblastn/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTN {
    take:
    assembly: Channel<Tuple<Map, Path>>
    query: Path

    main:
    TBLASTN_MODULE(assembly, query)
    CSVTK_CONCAT(gather(TBLASTN_MODULE.out, 'tblastn', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = TBLASTN_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
