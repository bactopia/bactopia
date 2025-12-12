/**
 * Determine the agr, spa and SCCmec types for _Staphylococcus aureus_ genomes.
 *
 * This subworkflow performs comprehensive typing of *Staphylococcus aureus* genomes by
 * determining the agr locus type using [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE),
 * spa repeat type using [spaTyper](https://github.com/HCGB-IGTP/spaTyper), and SCCmec element
 * type using SCCmec typing. It combines results from multiple typing methods to provide
 * a complete characterization of *S. aureus* strains.
 *
 * @status stable
 * @keywords staphylococcus aureus, agr typing, spa typing, sccmec, strain characterization
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation, database-dependent
 * @citation agrvate, spatyper, sccmec
 *
 * @subworkflows agrvate, spatyper, sccmec
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input repeats
 * Optional spa repeats database for improved spa typing
 *
 * @input repeat_order
 * Optional spa repeat order file for improved spa typing
 *
 * @output results      Aggregated results channel containing all typing results from agr, spa, and SCCmec analysis
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { AGRVATE      } from '../agrvate/main'
include { SPATYPER     } from '../spatyper/main'
include { SCCMEC       } from '../sccmec/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow STAPHTYPER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    repeats: Path?
    repeat_order: Path?

    main:
    // agrvate - agr locus type and agr operon variants
    AGRVATE(assembly)

    // spatyper - spa typing
    SPATYPER(assembly, repeats, repeat_order)

    // sccmec - SCCmec type based on targets and full cassettes
    SCCMEC(assembly)

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
