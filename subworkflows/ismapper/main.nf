//
// ismapper - Identify insertion sites positions in bacterial genomes
//
nextflow.preview.types = true

include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'

workflow ISMAPPER {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_reference  // channel: reference genome file
    ch_insertions // channel: insertion sequences file

    main:
    ISMAPPER_MODULE(ch_reads, ch_reference, ch_insertions)

    emit:
    results = ISMAPPER_MODULE.out.supplemental
    logs = ISMAPPER_MODULE.out.logs
    nf_logs = ISMAPPER_MODULE.out.nf_begin.mix(
        ISMAPPER_MODULE.out.nf_err,
        ISMAPPER_MODULE.out.nf_log,
        ISMAPPER_MODULE.out.nf_out,
        ISMAPPER_MODULE.out.nf_run,
        ISMAPPER_MODULE.out.nf_sh,
        ISMAPPER_MODULE.out.nf_trace
    )
    versions = ISMAPPER_MODULE.out.versions
}
