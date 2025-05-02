//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
include { ABRICATE_RUN } from '../../modules/abricate/run/main' //addParams( options: options )
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main' //addParams( options: [logs_subdir: 'abricate-concat', process_name: params.merge_folder] )

workflow ABRICATE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_abricate = Channel.empty()

    ABRICATE_RUN(fasta)
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions.first())

    ABRICATE_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:'abricate'], report]}.set{ ch_merge_abricate }
    ABRICATE_SUMMARY(ch_merge_abricate)
    ch_merged_abricate = ch_merged_abricate.mix(ABRICATE_SUMMARY.out.report)
    ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)

    emit:
    tsv = ABRICATE_RUN.out.report
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
    merged_tsv = ch_merged_abricate
    versions = ch_versions // channel: [ versions.yml ]
}
