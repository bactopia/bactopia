//
// panaroo - Pipeline for pangenome investigations
//
include { PANAROO_RUN } from '../../modules/panaroo/run/main'

workflow PANAROO {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    PANAROO_RUN(gff)

    emit:
    // Individual outputs
    csv = PANAROO_RUN.out.csv
    aln = PANAROO_RUN.out.aln
    filtered_aln = PANAROO_RUN.out.filtered_aln
    panaroo_csv = PANAROO_RUN.out.panaroo_csv

    // Generic aggregate outputs
    results = PANAROO_RUN.out.supplemental.mix(
        PANAROO_RUN.out.csv,
        PANAROO_RUN.out.aln,
        PANAROO_RUN.out.filtered_aln,
        PANAROO_RUN.out.panaroo_csv
    )
    logs = PANAROO_RUN.out.logs
    nf_logs = PANAROO_RUN.out.nf_begin.mix(
        PANAROO_RUN.out.nf_err,
        PANAROO_RUN.out.nf_log,
        PANAROO_RUN.out.nf_out,
        PANAROO_RUN.out.nf_run,
        PANAROO_RUN.out.nf_sh,
        PANAROO_RUN.out.nf_trace
    )
    versions = PANAROO_RUN.out.versions
}
