//
// amrfinderplus - Identify antimicrobial resistance in genes or proteins
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus')
options.args = [
    params.report_common ? "--report_common" : "",
    params.report_all_equal ? "--report_all_equal" : "",
    params.organism ? "--organism ${params.organism}" : "",
    params.amrfinder_noplus ? "" : "--plus",
    "--ident_min ${params.ident_min}",
    "--coverage_min ${params.coverage_min}",
    "--translation_table ${params.translation_table}",
    "${params.amrfinder_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { AMRFINDERPLUS_RUN } from '../../../modules/nf-core/amrfinderplus/run/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'amrfinderplus-concat', process_name: params.merge_folder] )

workflow AMRFINDERPLUS {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ amrfinderplus_db ]

    main:
    ch_versions = Channel.empty()
    ch_amrfinder_db = Channel.empty()
    ch_merged_reports = Channel.empty()

    // Run AMRFinder
    if (params.amrfinder_db) {
        AMRFINDERPLUS_RUN(fasta, file(params.amrfinder_db))
    } else {
        AMRFINDERPLUS_RUN(fasta, db)
    }
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())

    // Merge results
    AMRFINDERPLUS_RUN.out.report.collect{meta, report -> report}.map{ report -> [[id:'amrfinderplus'], report]}.set{ ch_merge_report }
    CSVTK_CONCAT(ch_merge_report, 'tsv', 'tsv')
    ch_merged_reports = ch_merged_reports.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    report = AMRFINDERPLUS_RUN.out.report
    merged_tsv = ch_merged_reports
    mutation_report = AMRFINDERPLUS_RUN.out.mutation_report
    db = ch_amrfinder_db
    versions = ch_versions // channel: [ versions.yml ]
}
