//
// scoary - GWAS analysis using pangenome outputs
//
include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'

workflow SCOARY {
    take:
    csv // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    SCOARY_TRAITS = params.traits ? file(params.traits) : []
    
    SCOARY_MODULE(csv, SCOARY_TRAITS)
    ch_versions = ch_versions.mix(SCOARY_MODULE.out.versions)

    emit:
    csv = SCOARY_MODULE.out.csv
    logs = SCOARY_MODULE.out.logs
    nf_logs = SCOARY_MODULE.out.nf_begin.mix(
        SCOARY_MODULE.out.nf_err,
        SCOARY_MODULE.out.nf_log,
        SCOARY_MODULE.out.nf_out,
        SCOARY_MODULE.out.nf_run,
        SCOARY_MODULE.out.nf_sh,
        SCOARY_MODULE.out.nf_trace
    )
    versions = ch_versions
}
