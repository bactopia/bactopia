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
    ch_nf_logs = Channel.empty()

    PANAROO_RUN(gff)
    ch_versions = ch_versions.mix(PANAROO_RUN.out.versions)
    ch_logs = ch_logs.mix(PANAROO_RUN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PANAROO_RUN.out.nf_logs)

    emit:
    aln = PANAROO_RUN.out.aln
    csv = PANAROO_RUN.out.csv
    panaroo_csv = PANAROO_RUN.out.panaroo_csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
