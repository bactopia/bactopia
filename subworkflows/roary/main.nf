//
// roary - Rapid large-scale prokaryote pangenome analysis
//
include { ROARY as ROARY_MODULE } from '../../modules/roary/main'

workflow ROARY {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ROARY_MODULE(gff)
    ch_versions = ch_versions.mix(ROARY_MODULE.out.versions)

    emit:
    aln = ROARY_MODULE.out.aln
    csv = ROARY_MODULE.out.csv
    logs = ROARY_MODULE.out.logs
    nf_logs = ROARY_MODULE.out.nf_begin.mix(
        ROARY_MODULE.out.nf_err,
        ROARY_MODULE.out.nf_log,
        ROARY_MODULE.out.nf_out,
        ROARY_MODULE.out.nf_run,
        ROARY_MODULE.out.nf_sh,
        ROARY_MODULE.out.nf_trace
    )
    versions = ch_versions
}
