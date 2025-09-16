//
// ismapper - Identify insertion sites positions in bacterial genomes
//
include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'

workflow ISMAPPER {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_reference  // channel: reference genome file
    ch_insertions // channel: insertion sequences file

    main:
    ch_versions = Channel.empty()

    ISMAPPER_MODULE(ch_reads, ch_reference, ch_insertions)
    ch_versions = ch_versions.mix(ISMAPPER_MODULE.out.versions)

    emit:
    results = ISMAPPER_MODULE.out.results
    logs = ISMAPPER_MODULE.out.logs
    nf_logs = ISMAPPER_MODULE.out.nf_begin.mix(
        ISMAPPER_MODULE.out.nf_err,
        ISMAPPER_MODULE.out.nf_log,
        ISMAPPER_MODULE.out.nf_out,
        ISMAPPER_MODULE.out.nf_run,
        ISMAPPER_MODULE.out.nf_sh,
        ISMAPPER_MODULE.out.nf_trace
    )
    versions = ch_versions
}
