/**
 * Calculate Mash distances between sequences and a reference.
 *
 * This subworkflow uses [Mash](https://github.com/marbl/Mash) to calculate MinHash-based
 * distances between query sequences and a reference sequence. It creates Mash sketches
 * of the input sequences and computes distance values, then aggregates all distance
 * calculations into a single consolidated report.
 *
 * @status stable
 * @keywords mash, distance, minhash, comparison, reference
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation mash
 *
 * @modules csvtk_concat, mash_dist
 *
 * @input tuple(meta, seqs)
 * - `meta`: Groovy Map containing sample information
 * - `seqs`: Sequences in FASTA format to compare against reference
 *
 * @input reference
 * Reference sequence in FASTA format for distance calculations
 *
 * @output dist            Per-sample TSV files containing Mash distance calculations
 * @output merged_dist     Consolidated TSV file containing Mash distances from all samples
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MASH_DIST    } from '../../modules/mash/dist/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow MASHDIST {
    take:
    seqs: Channel<Tuple<Map, Set<Path>>>
    reference: Path

    main:
    MASH_DIST(seqs, reference)
    CSVTK_CONCAT(gather(MASH_DIST.out.dist, 'mashdist', 'dist'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    dist: Channel<Tuple<Map, Set<Path>>> = MASH_DIST.out.dist
    merged_dist: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.dist,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MASH_DIST.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
