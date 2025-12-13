/**
 * Predict serotypes of Shigella and EIEC from assemblies.
 *
 * This subworkflow uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to predict
 * serotypes of *Shigella* and Enteroinvasive *E. coli* (EIEC) from assembled genomes.
 * It uses a cluster-informed approach to identify specific serotype markers and classify
 * isolates based on their antigenic profiles.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, typing, cluster analysis
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation shigeifinder
 *
 * @modules shigeifinder, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing Shigella/EIEC serotype predictions
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../modules/shigeifinder/main'
include { CSVTK_CONCAT                        } from '../../modules/csvtk/concat/main'
include { flattenPaths                        } from 'plugin/nf-bactopia'
include { gather                              } from 'plugin/nf-bactopia'

workflow SHIGEIFINDER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    SHIGEIFINDER_MODULE(assembly)
    CSVTK_CONCAT(gather(SHIGEIFINDER_MODULE.out.tsv, 'shigeifinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SHIGEIFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
