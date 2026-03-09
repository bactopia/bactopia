/**
 * Search protein sequences against protein database.
 *
 * This subworkflow uses [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastp)
 * from the NCBI BLAST+ suite to search protein sequences against a protein database.
 * It processes each assembly individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords blast, protein, alignment, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation blast
 *
 * @modules blastp, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format that will be translated and searched
 *
 * @input query
 * Path to protein database for searching against translated sequences
 *
 * @output sample_outputs  Record outputs from BLASTP module per sample
 * @output run_outputs   Record with merged CSV and associated logs/versions
 */
nextflow.preview.types = true

include { BLAST_BLASTP as BLASTP_MODULE } from '../../modules/blast/blastp/main'
include { CSVTK_CONCAT                 } from '../../modules/csvtk/concat/main'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTP {
    take:
    assembly: Channel<Tuple<Map, Path>>
    query: Path

    main:
    BLASTP_MODULE(assembly, query)
    CSVTK_CONCAT(gather(BLASTP_MODULE.out, 'blastp', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = BLASTP_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
