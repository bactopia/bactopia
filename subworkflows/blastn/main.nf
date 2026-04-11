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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation blast
 *
 * @modules blast_blastn, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
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
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow BLASTN {
    take:
    blastdb: Channel<Record>
    query: Value<Path>

    main:
    ch_blastn = BLASTN_MODULE(blastdb, query)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_blastn, 'tsv', [name: 'blastn']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_blastn
    run_outputs = ch_csvtk_concat
}
