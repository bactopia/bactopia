//
// mykrobe - Antimicrobial resistance detection for specific species
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mykrobe')
options.is_module = params.wf == 'mykrobe' ? true : false
options.args = [
    params.report_all_calls ? "--report_all_calls" : "",
    "--kmer ${params.kmer}",
    "--min_depth ${params.min_depth}",
    "--model ${params.model}",
    "${params.mykrobe_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MYKROBE_PREDICT }  from '../../../modules/nf-core/modules/mykrobe/predict/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'mykrobe'] )
}

workflow MYKROBE {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_mykrobe = Channel.empty()

    MYKROBE_PREDICT(reads, params.mykrobe_species)
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())

    if (params.is_subworkflow) {
        MYKROBE_PREDICT.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'mykrobe'], csv]}.set{ ch_merge_mykrobe }
        CSVTK_CONCAT(ch_merge_mykrobe, 'csv', 'csv')
        ch_merged_mykrobe = ch_merged_mykrobe.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    csv = MYKROBE_PREDICT.out.csv
    merged_csv = ch_merged_mykrobe
    json = MYKROBE_PREDICT.out.json
    versions = ch_versions // channel: [ versions.yml ]
}
