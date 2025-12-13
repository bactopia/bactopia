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
 * @output tsv         Per-sample TSV files containing BLASTX alignment results
 * @output merged_tsv  Consolidated TSV file containing BLASTX results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BLAST_BLASTX as BLASTX_MODULE } from '../../modules/blast/blastx/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTX {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    BLASTX_MODULE(assembly, query)
    CSVTK_CONCAT(gather(BLASTX_MODULE.out.tsv, 'blastx'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = BLASTX_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTX_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTX_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTX_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTX_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
