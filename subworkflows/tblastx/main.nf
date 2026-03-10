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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation blast
 *
 * @modules tblastx, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format that will be translated and searched
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
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTX {
    take:
    assembly: Channel<Record>
    query: Path

    main:
    TBLASTX_MODULE(assembly, query)
    CSVTK_CONCAT(gather(TBLASTX_MODULE.out, 'tblastx', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = TBLASTX_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
