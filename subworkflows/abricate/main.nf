//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
include { ABRICATE_RUN } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'

workflow ABRICATE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ABRICATE_RUN(fasta)

    // Merge results
    ABRICATE_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:'abricate'], report]}.set{ ch_merge_abricate }
    ABRICATE_SUMMARY(ch_merge_abricate)

    emit:
    // Individual outputs
    tsv = ABRICATE_RUN.out.report
    merged_tsv = ABRICATE_SUMMARY.out.report

    // Generic aggregate outputs
    results = ABRICATE_RUN.out.report.mix(
        ABRICATE_SUMMARY.out.report
    )
    logs = ABRICATE_RUN.out.logs.mix(
        ABRICATE_SUMMARY.out.logs
    )
    nf_logs = ABRICATE_RUN.out.nf_begin.mix(
        ABRICATE_RUN.out.nf_err,
        ABRICATE_RUN.out.nf_log,
        ABRICATE_RUN.out.nf_out,
        ABRICATE_RUN.out.nf_run,
        ABRICATE_RUN.out.nf_sh,
        ABRICATE_RUN.out.nf_trace,
        ABRICATE_SUMMARY.out.nf_begin,
        ABRICATE_SUMMARY.out.nf_err,
        ABRICATE_SUMMARY.out.nf_log,
        ABRICATE_SUMMARY.out.nf_out,
        ABRICATE_SUMMARY.out.nf_run,
        ABRICATE_SUMMARY.out.nf_sh,
        ABRICATE_SUMMARY.out.nf_trace
    )
    versions = ABRICATE_RUN.out.versions.mix(
        ABRICATE_SUMMARY.out.versions
    )
}
