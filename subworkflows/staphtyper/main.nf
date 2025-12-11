/**
 * Determine the agr, spa and SCCmec types for _Staphylococcus aureus_ genomes.
 *
 * This subworkflow orchestrates the execution of main.nf analysis components.
 *
 * @status stable
 * @keywords main.nf, subworkflow, analysis
 * @tags complexity:simple input-type:multiple output-type:multiple features:components
 * @citation main.nf
 *
 *
 * @input fasta
 * Channel containing tuples with metadata and file paths
 *
 * @input repeats
 * Input channel
 *
 * @input repeat_order
 * Input channel
 *
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { AGRVATE      } from '../agrvate/main'
include { SPATYPER     } from '../spatyper/main'
include { SCCMEC       } from '../sccmec/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow STAPHTYPER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    repeats: Path?
    repeat_order: Path?

    main:
    // agrvate - agr locus type and agr operon variants
    AGRVATE(fasta)

    // spatyper - spa typing
    SPATYPER(fasta, repeats, repeat_order)

    // sccmec - SCCmec type based on targets and full cassettes





    SCCMEC(fasta)

    emit:
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE.out.results,
        SPATYPER.out.results,
        SCCMEC.out.results
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE.out.logs,
        SPATYPER.out.logs,
        SCCMEC.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE.out.nf_logs,
        SPATYPER.out.nf_logs,
        SCCMEC.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE.out.versions,
        SPATYPER.out.versions,
        SCCMEC.out.versions
    ])
}
