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
 * @output tsv         Per-sample TSV files containing TBLASTX alignment results
 * @output merged_tsv  Consolidated TSV file containing TBLASTX results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BLAST_TBLASTX as TBLASTX_MODULE } from '../../modules/blast/tblastx/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTX {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    TBLASTX_MODULE(assembly, query)
    CSVTK_CONCAT(gather(TBLASTX_MODULE.out.tsv, 'tblastx'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = TBLASTX_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTX_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTX_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTX_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTX_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
