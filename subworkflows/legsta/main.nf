//
// legsta - Typing of Legionella pneumophila assemblies
//
include { LEGSTA as LEGSTA_MODULE } from '../../modules/legsta/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LEGSTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    LEGSTA_MODULE(fasta)
    ch_versions = ch_versions.mix(LEGSTA_MODULE.out.versions)
    ch_logs = ch_logs.mix(LEGSTA_MODULE.out.logs)

    // Merge results
    LEGSTA_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'legsta'], tsv]}.set{ ch_merge_legsta }
    CSVTK_CONCAT(ch_merge_legsta, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = LEGSTA_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        LEGSTA_MODULE.out.nf_begin,
        LEGSTA_MODULE.out.nf_err,
        LEGSTA_MODULE.out.nf_log,
        LEGSTA_MODULE.out.nf_out,
        LEGSTA_MODULE.out.nf_run,
        LEGSTA_MODULE.out.nf_sh,
        LEGSTA_MODULE.out.nf_trace
    )
    versions = ch_versions
}
