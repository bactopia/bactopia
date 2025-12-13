/**
 * Predict serogroups of Pseudomonas aeruginosa from assemblies.
 *
 * This subworkflow uses [Pasty](https://github.com/rpetit3/pasty) to perform in silico
 * serogrouping of *Pseudomonas aeruginosa* isolates from assembled genomes. It identifies
 * O-antigen biosynthesis genes to classify isolates into their known serogroups using
 * BLAST-based homology searches.
 *
 * @status stable
 * @keywords pseudomonas aeruginosa, serogroup, typing, o-antigen, prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation pasty
 *
 * @modules pasty, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing P. aeruginosa serogroup predictions
 * @output merged_tsv  Consolidated TSV file containing serogroup predictions from all samples
 * @output blast       Per-sample BLAST results showing hits to O-antigen genes
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow PASTY {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    PASTY_MODULE(assembly)
    CSVTK_CONCAT(gather(PASTY_MODULE.out.tsv, 'pasty'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = PASTY_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    blast: Channel<Tuple<Map, Set<Path>>> = PASTY_MODULE.out.blast

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.tsv,
        PASTY_MODULE.out.blast,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
