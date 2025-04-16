//
// abricate - Mass screening of contigs for antimicrobial and virulence genes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'abricate')
options.args = [
    "--db ${params.abricate_db}",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
options.subdir = params.abricate_db
include { ABRICATE_RUN } from '../../../modules/nf-core/abricate/run/main' addParams( options: options )
include { ABRICATE_SUMMARY } from '../../../modules/nf-core/abricate/summary/main' addParams( options: [logs_subdir: 'abricate-concat', process_name: params.merge_folder] )

workflow ABRICATE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_abricate = Channel.empty()

    ABRICATE_RUN(fasta)
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions.first())

    ABRICATE_RUN.out.report.collect{meta, report -> report}.map{ report -> [[id:'abricate'], report]}.set{ ch_merge_abricate }
    ABRICATE_SUMMARY(ch_merge_abricate)
    ch_merged_abricate = ch_merged_abricate.mix(ABRICATE_SUMMARY.out.report)
    ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)

    emit:
    tsv = ABRICATE_RUN.out.report
    merged_tsv = ch_merged_abricate
    versions = ch_versions // channel: [ versions.yml ]
}
