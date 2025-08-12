//
// panaroo - Pipeline for pangenome investigations
//
include { PANAROO_RUN } from '../../modules/panaroo/run/main'

workflow PANAROO {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    PANAROO_RUN(gff)
    ch_versions = ch_versions.mix(PANAROO_RUN.out.versions)
    ch_logs = ch_logs.mix(PANAROO_RUN.out.logs)

    emit:
    csv = PANAROO_RUN.out.csv
    aln = PANAROO_RUN.out.aln
    panaroo_csv = PANAROO_RUN.out.panaroo_csv
    logs = ch_logs
    nf_logs = PANAROO_RUN.out.nf_begin.mix(
        PANAROO_RUN.out.nf_err,
        PANAROO_RUN.out.nf_log,
        PANAROO_RUN.out.nf_out,
        PANAROO_RUN.out.nf_run,
        PANAROO_RUN.out.nf_sh,
        PANAROO_RUN.out.nf_trace
    )
    versions = ch_versions
}
