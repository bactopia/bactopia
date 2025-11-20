//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
nextflow.preview.types = true

include { ABRICATE_RUN } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'

workflow ABRICATE {
    take:
    fasta: Channel<Tuple<Map,Path>>

    main:
    ABRICATE_RUN(fasta)

    // Merge results
    ABRICATE_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:'abricate'], report]}.set{ ch_merge_abricate }
    ABRICATE_SUMMARY(ch_merge_abricate)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map,Path>> = ABRICATE_RUN.out.report
    merged_tsv: Channel<Tuple<Map,Path>> = channel.of(ABRICATE_SUMMARY.out.report)

    // Generic aggregate outputs
    results: Channel<Tuple<Map,Path>> = ABRICATE_RUN.out.report.mix(
        channel.of(ABRICATE_SUMMARY.out.report)
    )
    logs: Channel<Tuple<Map,Path>> = ABRICATE_RUN.out.logs.mix(
        channel.of(ABRICATE_SUMMARY.out.logs)
    )
    nf_logs: Channel<Tuple<Map,Path>> = ABRICATE_RUN.out.nf_begin.mix(
        ABRICATE_RUN.out.nf_err,
        ABRICATE_RUN.out.nf_log,
        ABRICATE_RUN.out.nf_out,
        ABRICATE_RUN.out.nf_run,
        ABRICATE_RUN.out.nf_sh,
        ABRICATE_RUN.out.nf_trace,
        channel.of(ABRICATE_SUMMARY.out.nf_begin),
        channel.of(ABRICATE_SUMMARY.out.nf_err),
        channel.of(ABRICATE_SUMMARY.out.nf_log),
        channel.of(ABRICATE_SUMMARY.out.nf_out),
        channel.of(ABRICATE_SUMMARY.out.nf_run),
        channel.of(ABRICATE_SUMMARY.out.nf_sh),
        channel.of(ABRICATE_SUMMARY.out.nf_trace)
    )
    versions: Channel<Tuple<Map,Path>> = ABRICATE_RUN.out.versions.mix(
        channel.of(ABRICATE_SUMMARY.out.versions)
    )
}
