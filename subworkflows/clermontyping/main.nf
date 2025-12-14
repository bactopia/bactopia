/**
 * Predict phylogroups of Escherichia coli from genome assemblies.
 *
 * This subworkflow uses [ClermontTyping](https://github.com/happykhan/ClermonTyping) to determine
 * the phylogenetic groups of *Escherichia coli* strains from assembled genomes. It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords escherichia coli, phylogroup, typing, clermont
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation clermontyping
 *
 * @modules csvtk_concat, clermontyping
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv            Per-sample TSV files containing ClermontTyping results
 * @output merged_tsv     Consolidated TSV file containing ClermontTyping from all samples
 * @output supplemental   Additional supplemental output files from ClermontTyping
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../modules/clermontyping/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { flattenPaths                          } from 'plugin/nf-bactopia'
include { gather                                } from 'plugin/nf-bactopia'

workflow CLERMONTYPING {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    CLERMONTYPING_MODULE(assembly)
    CSVTK_CONCAT(gather(CLERMONTYPING_MODULE.out.tsv, 'clermontyping'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = CLERMONTYPING_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        CLERMONTYPING_MODULE.out.tsv,
        CLERMONTYPING_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CLERMONTYPING_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CLERMONTYPING_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        CLERMONTYPING_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
