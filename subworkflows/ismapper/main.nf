/**
 * Identify insertion sites positions in bacterial genomes.
 *
 * This subworkflow orchestrates the execution of main.nf analysis components.
 *
 * @status stable
 * @keywords main.nf, subworkflow, analysis
 * @tags complexity:simple input-type:multiple output-type:multiple
 * @citation main.nf
 *
 * @modules ismapper as ismapper_module
 *
 * @input ch_reads
 * Channel containing tuples with metadata and file paths
 *
 * @input ch_reference
 * Input channel
 *
 * @input ch_insertions
 * Input channel
 *
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */


nextflow.preview.types = true

include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow ISMAPPER {
    take:
    ch_reads: Channel<Tuple<Map, Set<Path>>>
    ch_reference: Path
    ch_insertions: Path

    main:
    ISMAPPER_MODULE(ch_reads, ch_reference, ch_insertions)

    emit:
    results: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.supplemental
    logs: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.versions
}
