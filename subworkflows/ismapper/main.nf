//
// ismapper - Identify insertion sites positions in bacterial genomes
//
include { ISMAPPER } from '../../modules/ismapper/main'

workflow ISMAPPER_WORKFLOW {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_reference  // channel: reference genome file
    ch_insertions // channel: insertion sequences file

    main:
    ch_versions = Channel.empty()

    ISMAPPER(ch_reads, ch_reference, ch_insertions)
    ch_versions = ch_versions.mix(ISMAPPER.out.versions.first())

    emit:
    results  = ISMAPPER.out.results
    logs     = ISMAPPER.out.logs
    nf_logs  = ISMAPPER.out.nf_begin.mix(ISMAPPER.out.nf_err, ISMAPPER.out.nf_log, ISMAPPER.out.nf_out, ISMAPPER.out.nf_run, ISMAPPER.out.nf_sh, ISMAPPER.out.nf_trace)
    versions = ch_versions
}
