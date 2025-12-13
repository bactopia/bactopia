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
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format to be searched against the database
 *
 * @input query
 * FASTA file containing nucleotide query sequences to search for in the assemblies
 *
 * @output tsv         Per-sample tab-delimited summaries of alignments
 * @output merged_tsv  Consolidated tab-delimited file containing alignments from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BLAST_BLASTN as BLASTN_MODULE } from '../../modules/blast/blastn/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTN {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    BLASTN_MODULE(assembly, query)
    CSVTK_CONCAT(gather(BLASTN_MODULE.out.tsv, 'blastn'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = BLASTN_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
