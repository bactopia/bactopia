//
// staphtyper - Determine the agr, spa and SCCMEC types for S. aureus assemblies
//
nextflow.preview.types = true

include { AGRVATE      } from '../agrvate/main'
include { SPATYPER     } from '../spatyper/main'
include { SCCMEC       } from '../sccmec/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow STAPHTYPER {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ assemblies ] ]
    repeats: Channel<Tuple<Map, Path>>
    repeat_order: Channel<Tuple<Map, Path>>

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
