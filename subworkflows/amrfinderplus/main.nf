//
// amrfinderplus - Identify antimicrobial resistance in genes or proteins
//
nextflow.preview.types = true

include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow AMRFINDERPLUS {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ amrfinderplus_db ]

    main:
    // Run AMRFinder
    AMRFINDERPLUS_RUN(fasta, db)

    // Merge results
    ch_merge_report = AMRFINDERPLUS_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:'amrfinderplus'], report]}
    CSVTK_CONCAT(ch_merge_report, 'tsv', 'tsv')

    emit:
    // Individual outputs
    report = AMRFINDERPLUS_RUN.out.report
    merged_tsv = CSVTK_CONCAT.out.csv
    mutation_report = AMRFINDERPLUS_RUN.out.mutation_report

    // Generic aggregate outputs
    results = AMRFINDERPLUS_RUN.out.report.mix(
        AMRFINDERPLUS_RUN.out.mutation_report,
        CSVTK_CONCAT.out.csv
    )
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
    versions = AMRFINDERPLUS_RUN.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
