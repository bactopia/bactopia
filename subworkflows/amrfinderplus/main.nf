//
// amrfinderplus - Identify antimicrobial resistance in genes or proteins
//
include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

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
    AMRFINDERPLUS_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:'amrfinderplus'], report]}.set{ ch_merge_report }
    CSVTK_CONCAT(ch_merge_report, 'tsv', 'tsv')
    ch_merged_reports = ch_merged_reports.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    report = AMRFINDERPLUS_RUN.out.report
    merged_tsv = ch_merged_reports
    mutation_report = AMRFINDERPLUS_RUN.out.mutation_report
    db = ch_amrfinder_db
    logs = AMRFINDERPLUS_RUN.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = AMRFINDERPLUS_RUN.out.nf_begin.mix(
        AMRFINDERPLUS_RUN.out.nf_err,
        AMRFINDERPLUS_RUN.out.nf_log,
        AMRFINDERPLUS_RUN.out.nf_out,
        AMRFINDERPLUS_RUN.out.nf_run,
        AMRFINDERPLUS_RUN.out.nf_sh,
        AMRFINDERPLUS_RUN.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
