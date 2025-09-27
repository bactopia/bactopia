//
// checkm2 - Assess the assembly quality of your samples
//
include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CHECKM2 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    database // channel: [ database ]
    download_checkm2 // boolean

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_merged_checkm2 = Channel.empty()

    if (download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(fasta, CHECKM2_DOWNLOAD.out.db)
    } else {
        CHECKM2_PREDICT(fasta, database)
    }
    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)
    ch_logs = ch_logs.mix(CHECKM2_PREDICT.out.logs)

    // Merge results
    CHECKM2_PREDICT.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'checkm2'], tsv]}.set{ ch_merge_checkm2 }
    CSVTK_CONCAT(ch_merge_checkm2, 'tsv', 'tsv')
    ch_merged_checkm2 = ch_merged_checkm2.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    report = CHECKM2_PREDICT.out.tsv
    merged_reports = ch_merged_checkm2
    results = CHECKM2_PREDICT.out.results
    logs = ch_logs
    nf_logs = CHECKM2_PREDICT.out.nf_begin.mix(
        CHECKM2_PREDICT.out.nf_err,
        CHECKM2_PREDICT.out.nf_log,
        CHECKM2_PREDICT.out.nf_out,
        CHECKM2_PREDICT.out.nf_run,
        CHECKM2_PREDICT.out.nf_sh,
        CHECKM2_PREDICT.out.nf_trace,
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
