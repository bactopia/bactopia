//
// checkm2 - Assess the assembly quality of your samples
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'checkm2')
options.is_module = params.wf == 'checkm2' ? true : false
options.args = [
    params.checkm2_general ? "--general" : "",
    params.checkm2_specific ? "--specific" : "",
    params.checkm2_allmodels ? "--allmodels" : "",
    params.checkm2_genes ? "--genes" : "",
    "${params.checkm2_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE = params.checkm2_db ? file(params.checkm2_db) : []


include { CHECKM2_DOWNLOAD } from '../../../modules/nf-core/checkm2/download/main' addParams( )
include { CHECKM2_PREDICT } from '../../../modules/nf-core/checkm2/predict/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'checkm2-concat', process_name: params.merge_folder] )

workflow CHECKM2 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_checkm2 = Channel.empty()

    if (params.download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(fasta, CHECKM2_DOWNLOAD.out.db)
        ch_versions = ch_versions.mix(CHECKM2_DOWNLOAD.out.versions)
    } else {
        CHECKM2_PREDICT(fasta, DATABASE)
    }

    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())

    CHECKM2_PREDICT.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'checkm2'], tsv]}.set{ ch_merge_checkm2 }
    CSVTK_CONCAT(ch_merge_checkm2, 'tsv', 'tsv')
    ch_merged_checkm2 = ch_merged_checkm2.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    results = CHECKM2_PREDICT.out.results
    report = CHECKM2_PREDICT.out.tsv
    merged_reports = ch_merged_checkm2
    versions = ch_versions // channel: [ versions.yml ]
}
