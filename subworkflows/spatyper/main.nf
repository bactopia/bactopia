//
// SpaTyper - Computational method for finding spa types in Staphylococcus aureus
//
include { SPATYPER as SPATYPER_MODULE } from '../../modules/spatyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SPATYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    REPEATS = params.repeats ? file(params.repeats, checkIfExists: true) : []
    REPEAT_ORDER = params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []

    SPATYPER_MODULE(fasta, REPEATS, REPEAT_ORDER)
    ch_versions = ch_versions.mix(SPATYPER_MODULE.out.versions.first())

    // Merge results
    SPATYPER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'spatyper'], tsv]}.set{ ch_merge_spatyper }
    CSVTK_CONCAT(ch_merge_spatyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SPATYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = SPATYPER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SPATYPER_MODULE.out.nf_begin.mix(
        SPATYPER_MODULE.out.nf_err,
        SPATYPER_MODULE.out.nf_log,
        SPATYPER_MODULE.out.nf_out,
        SPATYPER_MODULE.out.nf_run,
        SPATYPER_MODULE.out.nf_sh,
        SPATYPER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions
}
