//
// shigapass - Predict Shigella serotypes and differentiate Shigella, EIEC and non-Shigella/EIEC
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigapass')

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../../modules/nf-core/shigapass/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'shigapass-summary-concat', process_name: params.merge_folder] )

workflow SHIGAPASS {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SHIGAPASS_MODULE(fasta)
    ch_versions = ch_versions.mix(SHIGAPASS_MODULE.out.versions.first())

    // Merge results
    SHIGAPASS_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'shigapass'], tsv]}.set{ ch_merge_shigapass }
    CSVTK_CONCAT(ch_merge_shigapass, 'tsv', 'tsv')

    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SHIGAPASS_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
