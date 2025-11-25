//
// checkm2 - Assess the assembly quality of your samples
//
nextflow.preview.types = true

include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CHECKM2 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    database // channel: [ database ]
    download_checkm2 // boolean

    main:
    if (download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(fasta, CHECKM2_DOWNLOAD.out.db)
    } else {
        CHECKM2_PREDICT(fasta, database)
    }

    // Merge results
    ch_merge_checkm2 = CHECKM2_PREDICT.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'checkm2'], tsv]}
    CSVTK_CONCAT(ch_merge_checkm2, 'tsv', 'tsv')

    emit:
    // Individual outputs
    report = CHECKM2_PREDICT.out.tsv
    merged_reports = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = CHECKM2_PREDICT.out.tsv.mix(
        CHECKM2_PREDICT.out.supplemental,
        CSVTK_CONCAT.out.csv
    )
    logs = CHECKM2_PREDICT.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = CHECKM2_PREDICT.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
