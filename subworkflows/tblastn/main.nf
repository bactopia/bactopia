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
 * @output tsv         Per-sample TSV files containing TBLASTN alignment results
 * @output merged_tsv  Consolidated TSV file containing TBLASTN results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../modules/blast/tblastn/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTN {
    take:
    assembly: Channel<Tuple<Map, Path>>
    query: Path

    main:
    TBLASTN_MODULE(assembly, query)
    CSVTK_CONCAT(gather(TBLASTN_MODULE.out.tsv, 'tblastn'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = TBLASTN_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
