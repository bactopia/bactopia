//
// genotyphi - Salmonella Typhi genotyping with Mykrobe outputs
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'genotyphi')
mykrobe_args = [
    params.report_all_calls ? "--report_all_calls" : "",
    "--kmer ${params.kmer}",
    "--min_depth ${params.min_depth}",
    "--model ${params.model}",
    "${params.mykrobe_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MYKROBE_PREDICT } from '../../../modules/nf-core/mykrobe/predict/main' addParams( options: options + [ args: "${mykrobe_args}", logs_subdir: 'mykrobe'] )
include { GENOTYPHI_PARSE  } from '../../../modules/nf-core/genotyphi/parse/main' addParams( options: options + [ args: "${mykrobe_args}", logs_subdir: 'genotyphi'] )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'genotyphi-concat', process_name: params.merge_folder] )

workflow GENOTYPHI {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    MYKROBE_PREDICT(reads, "typhi")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())

    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    ch_versions = ch_versions.mix(GENOTYPHI_PARSE.out.versions.first())

    // Merge results
    GENOTYPHI_PARSE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'genotyphi'], tsv]}.set{ ch_merge_genotyphi }
    CSVTK_CONCAT(ch_merge_genotyphi, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = GENOTYPHI_PARSE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    csv = MYKROBE_PREDICT.out.csv
    json = MYKROBE_PREDICT.out.json
    versions = ch_versions // channel: [ versions.yml ]
}
