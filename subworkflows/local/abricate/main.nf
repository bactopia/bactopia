//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'abricate')
options.is_module = params.wf == 'abricate' ? true : false
options.args = [
    "--db ${params.abricate_db}",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ABRICATE_RUN } from '../../../modules/nf-core/modules/abricate/run/main' addParams( options: options )
if (params.is_subworkflow) {
    include { ABRICATE_SUMMARY } from '../../../modules/nf-core/modules/abricate/summary/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'abricate'] )
}

workflow ABRICATE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_abricate = Channel.empty()

    ABRICATE_RUN(fasta)
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions.first())

    if (params.is_subworkflow) {
        ABRICATE_RUN.out.report.collect{meta, report -> report}.map{ report -> [[id:'abricate'], report]}.set{ ch_merge_abricate }
        ABRICATE_SUMMARY(ch_merge_abricate)
        ch_merged_abricate = ch_merged_abricate.mix(ABRICATE_SUMMARY.out.report)
        ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)
    }

    emit:
    report = ABRICATE_RUN.out.report
    merged_reports = ch_merged_abricate
    versions = ch_versions // channel: [ versions.yml ]
}
