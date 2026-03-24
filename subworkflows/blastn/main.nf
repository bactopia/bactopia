/**
 * Search a nucleotide database using nucleotide query sequences.
 *
 * This subworkflow uses [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to align
 * nucleotide query sequences against a nucleotide BLAST database. It processes each
 * sample individually and aggregates the alignment results into a single consolidated
 * report for all samples.
 *
 * @status stable
 * @keywords blast, alignment, nucleotide, search, fasta, database
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation blast
 *
 * @modules blastn, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format to be searched against the database
 *
 * @input query
 * FASTA file containing nucleotide query sequences to search for in the assemblies
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of alignments (standard BLAST outfmt 6)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BLAST_BLASTN as BLASTN_MODULE } from '../../modules/blast/blastn/main'
include { CSVTK_CONCAT                 } from '../../modules/csvtk/concat/main'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTN {
    take:
    assembly: Channel<Record>
    query: Path

    main:
    BLASTN_MODULE(assembly, query)
    CSVTK_CONCAT(gather(BLASTN_MODULE.out, 'tsv', [name: 'blastn']), 'tsv', 'tsv')

    emit:
    sample_outputs = BLASTN_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
