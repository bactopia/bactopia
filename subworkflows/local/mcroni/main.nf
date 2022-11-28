//
// mcroni - Sequence variation in mcr-1 genes (mobilized colistin resistance)
//

include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mcroni')
options.is_module = params.wf == 'mcroni' ? true : false
options.args = ""

include { MCRONI as MCRONI_MODULE } from '../../../modules/nf-core/mcroni/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'mcroni'] )
}

workflow MCRONI {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_mcroni = Channel.empty()

    MCRONI_MODULE(fasta)
    ch_versions = ch_versions.mix(MCRONI_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        MCRONI_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'mcroni'], tsv]}.set{ ch_merge_mcroni }
        CSVTK_CONCAT(ch_merge_mcroni, 'tsv', 'tsv')
        ch_merged_mcroni = ch_merged_mcroni.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = MCRONI_MODULE.out.tsv
    merged_tsv = ch_merged_mcroni
    fa = MCRONI_MODULE.out.fa
    versions = ch_versions
}
