//
// mcroni - Sequence variation in mcr-1 genes (mobilized colistin resistance)
//

include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mcroni')
options.args = ""

include { MCRONI as MCRONI_MODULE } from '../../../modules/nf-core/mcroni/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'mcroni'] )

workflow MCRONI {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MCRONI_MODULE(fasta)
    ch_versions = ch_versions.mix(MCRONI_MODULE.out.versions.first())

    // Merge results
    MCRONI_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'mcroni'], tsv]}.set{ ch_merge_mcroni }
    CSVTK_CONCAT(ch_merge_mcroni, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = MCRONI_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    fa = MCRONI_MODULE.out.fa
    versions = ch_versions
}
