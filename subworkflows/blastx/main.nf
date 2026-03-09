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
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format that will be translated and searched
 *
 * @input query
 * Path to BLAST protein database for searching translated sequences
 *
 * @output sample_outputs  Record outputs from BLASTX module per sample
 * @output run_outputs   Record with merged CSV and associated logs/versions
 */
nextflow.preview.types = true

include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTX {
    take:
    assembly: Channel<Tuple<Map, Path>>
    query: Path

    main:
    BLASTX_MODULE(assembly, query)
    CSVTK_CONCAT(gather(BLASTX_MODULE.out, 'blastx', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = BLASTX_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
