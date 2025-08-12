//
// tbprofiler - Detect resistance and lineages of Mycobacterium tuberculosis genomes
//
include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'

workflow TBPROFILER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    TBPROFILER_PROFILE(reads)
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE.out.versions)
    ch_logs = ch_logs.mix(TBPROFILER_PROFILE.out.logs)
    ch_versions = ch_versions.mix(TBPROFILER_COLLATE.out.versions)
    ch_logs = ch_logs.mix(TBPROFILER_COLLATE.out.logs)

    emit:
    csv = TBPROFILER_PROFILE.out.csv
    json = TBPROFILER_PROFILE.out.json
    txt = TBPROFILER_PROFILE.out.txt
    bam = TBPROFILER_PROFILE.out.bam
    vcf = TBPROFILER_PROFILE.out.vcf
    merged_csv = TBPROFILER_COLLATE.out.csv
    logs = ch_logs
    nf_logs = TBPROFILER_PROFILE.out.nf_begin.mix(
        TBPROFILER_PROFILE.out.nf_err,
        TBPROFILER_PROFILE.out.nf_log,
        TBPROFILER_PROFILE.out.nf_out,
        TBPROFILER_PROFILE.out.nf_run,
        TBPROFILER_PROFILE.out.nf_sh,
        TBPROFILER_PROFILE.out.nf_trace
    )
    versions = ch_versions
}
