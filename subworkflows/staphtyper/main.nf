//
// staphtyper - Determine the agr, spa and SCCMEC types for S. aureus assemblies
//
nextflow.preview.types = true

include { AGRVATE } from '../agrvate/main'
include { SPATYPER } from '../spatyper/main'
include { SCCMEC } from '../sccmec/main'

workflow STAPHTYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]
    repeats
    repeat_order

    main:
    // agrvate - agr locus type and agr operon variants
    AGRVATE(fasta)

    // spatyper - spa typing
    SPATYPER(fasta, repeats, repeat_order)

    // sccmec - SCCmec type based on targets and full cassettes
    SCCMEC(fasta)

    emit:
    results = AGRVATE.out.results.mix(
        SPATYPER.out.results,
        SCCMEC.out.results,
    )
    logs = AGRVATE.out.logs.mix(
        SPATYPER.out.logs,
        SCCMEC.out.logs,
    )
    nf_logs = AGRVATE.out.nf_logs.mix(
        SPATYPER.out.nf_logs,
        SCCMEC.out.nf_logs,
    )
    versions = AGRVATE.out.versions.mix(
        SPATYPER.out.versions,
        SCCMEC.out.versions,
    )
}
