//
// shigapass - Predict Shigella serotypes and differentiate Shigella, EIEC and non-Shigella/EIEC
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigapass')

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../../modules/nf-core/shigapass/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'shigapass-concat', process_name: params.merge_folder] )

workflow SHIGAPASS {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SHIGAPASS_MODULE(fasta)
    ch_versions = ch_versions.mix(SHIGAPASS_MODULE.out.versions.first())

    // Merge results
    SHIGAPASS_MODULE.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'shigapass'], csv]}.set{ ch_merge_shigapass }
    CSVTK_CONCAT(ch_merge_shigapass, 'csv', 'csv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    csv = SHIGAPASS_MODULE.out.csv
    merged_csv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
