//
// checkm2 - Assess the assembly quality of your samples
//
include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CHECKM2 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_merged_checkm2 = Channel.empty()

    DATABASE = params.checkm2_db ? file(params.checkm2_db) : []

    if (params.download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(fasta, CHECKM2_DOWNLOAD.out.db)
        ch_versions = ch_versions.mix(CHECKM2_DOWNLOAD.out.versions)
        ch_logs = ch_logs.mix(CHECKM2_DOWNLOAD.out.logs)
        ch_nf_logs = ch_nf_logs.mix(CHECKM2_DOWNLOAD.out.nf_logs)
    } else {
        CHECKM2_PREDICT(fasta, DATABASE)
    }

    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
    ch_logs = ch_logs.mix(CHECKM2_PREDICT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CHECKM2_PREDICT.out.nf_logs)

    CHECKM2_PREDICT.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'checkm2'], tsv] }.set{ ch_merge_checkm2 }
    CSVTK_CONCAT(ch_merge_checkm2, 'tsv', 'tsv')
    ch_merged_checkm2 = ch_merged_checkm2.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    results = CHECKM2_PREDICT.out.results
    report = CHECKM2_PREDICT.out.tsv
    merged_reports = ch_merged_checkm2
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
