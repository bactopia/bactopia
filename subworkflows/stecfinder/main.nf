/**
 * Identify and serotype Shiga toxin-producing E. coli (STEC) from assemblies.
 *
 * This subworkflow uses [STECFinder](https://github.com/LanLab/STECFinder) to identify
 * and serotype Shiga toxin-producing *E. coli* (STEC) strains using genomic cluster-specific
 * markers. It screens assemblies for virulence genes and serotype markers to classify
 * STEC isolates into their known serotypes.
 *
 * @status stable
 * @keywords escherichia coli, stec, serotype, virulence genes, shiga toxin
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation stecfinder
 *
 * @modules stecfinder, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing STEC identification and serotyping results
 * @output merged_tsv  Consolidated TSV file containing STEC results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { STECFINDER as STECFINDER_MODULE } from '../../modules/stecfinder/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow STECFINDER {
    take:
    seqs: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    STECFINDER_MODULE(seqs)
    CSVTK_CONCAT(gather(STECFINDER_MODULE.out.tsv, 'stecfinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = STECFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
