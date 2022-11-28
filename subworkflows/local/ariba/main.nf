//
// ariba - Gene identification through local assemblies
//

include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ariba')
options.is_module = params.wf == 'ariba' ? true : false
options.args = [
    params.spades_options ? "--spades_options '${params.spades_options}'" : "",
    params.ariba_no_clean ? "--noclean" : "",
    "--nucmer_min_id ${params.nucmer_min_id}",
    "--nucmer_min_len ${params.nucmer_min_len}",
    "--nucmer_breaklen ${params.nucmer_breaklen}",
    "--assembly_cov ${params.assembly_cov}",
    "--min_scaff_depth ${params.min_scaff_depth}",
    "--assembled_threshold ${params.assembled_threshold}",
    "--gene_nt_extend ${params.gene_nt_extend}",
    "--unique_threshold ${params.unique_threshold}",
    "--force",
    "--verbose",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ARIBA_GETREF } from '../../../modules/nf-core/ariba/getref/main' addParams( )
include { ARIBA_RUN } from '../../../modules/nf-core/ariba/run/main' addParams( options: options )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT as CSVTK_CONCAT_REPORT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [args: '-C "$" --lazy-quotes', publish_to_base: true, logs_subdir: 'ariba-report'] )
    include { CSVTK_CONCAT as CSVTK_CONCAT_SUMMARY } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [args: '--lazy-quotes', publish_to_base: true, logs_subdir: 'ariba-summary'] )
}

workflow ARIBA {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_report = Channel.empty()
    ch_merged_summary = Channel.empty()

    if (file("${params.ariba_dir}/${params.ariba_db}.tar.gz").exists()) {
        // Run Ariba
        ARIBA_RUN(reads, file("${params.ariba_dir}/${params.ariba_db}.tar.gz"))
        ch_versions = ch_versions.mix(ARIBA_RUN.out.versions.first())
    } else {
        // Build database and run Ariba
        ARIBA_GETREF(params.ariba_db)
        ARIBA_RUN(reads, ARIBA_GETREF.out.db)
        ch_versions = ch_versions.mix(ARIBA_GETREF.out.versions)
        ch_versions = ch_versions.mix(ARIBA_RUN.out.versions.first())
    }

    if (params.is_subworkflow) {
        ARIBA_RUN.out.report.collect{meta, report -> report}.map{ report -> [[id:'ariba-report'], report]}.set{ ch_merge_report }
        CSVTK_CONCAT_REPORT(ch_merge_report, 'tsv', 'tsv')
        ch_merged_report = ch_merged_report.mix(CSVTK_CONCAT_REPORT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT_REPORT.out.versions)

        ARIBA_RUN.out.summary.collect{meta, summary -> summary}.map{ summary -> [[id:'ariba-summary'], summary]}.set{ ch_merge_summary }
        CSVTK_CONCAT_SUMMARY(ch_merge_summary, 'csv', 'csv')
        ch_merged_summary = ch_merged_summary.mix(CSVTK_CONCAT_SUMMARY.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT_SUMMARY.out.versions)
    }

    emit:
    results = ARIBA_RUN.out.results
    merged_report = ch_merged_report
    merged_summary = ch_merged_summary
    versions = ch_versions // channel: [ versions.yml ]
}
