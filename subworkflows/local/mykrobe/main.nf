//
// mykrobe - Antimicrobial resistance detection for specific species
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mykrobe')
options.args = [
    params.report_all_calls ? "--report_all_calls" : "",
    "--kmer ${params.kmer}",
    "--min_depth ${params.min_depth}",
    "--model ${params.model}",
    "${params.mykrobe_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MYKROBE_PREDICT }  from '../../../modules/nf-core/mykrobe/predict/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'mykrobe-concat', process_name: params.merge_folder] )

workflow MYKROBE {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    MYKROBE_PREDICT(reads, params.mykrobe_species)
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())

    // Merge results
    MYKROBE_PREDICT.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'mykrobe'], csv]}.set{ ch_merge_mykrobe }
    CSVTK_CONCAT(ch_merge_mykrobe, 'csv', 'csv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    csv = MYKROBE_PREDICT.out.csv
    merged_csv = CSVTK_CONCAT.out.csv
    json = MYKROBE_PREDICT.out.json
    versions = ch_versions // channel: [ versions.yml ]
}
