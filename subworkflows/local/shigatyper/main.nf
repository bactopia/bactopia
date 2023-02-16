//
// shigatyper - Shigella serotype from Illumina or Oxford Nanopore reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigatyper')

include { SHIGATYPER as SHIGATYPER_MODULE } from '../../../modules/nf-core/shigatyper/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'shigatyper'] )

workflow SHIGATYPER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    SHIGATYPER_MODULE(reads)
    ch_versions = ch_versions.mix(SHIGATYPER_MODULE.out.versions.first())

    // Merge results
    SHIGATYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'shigatyper'], tsv]}.set{ ch_merge_shigatyper }
    CSVTK_CONCAT(ch_merge_shigatyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SHIGATYPER_MODULE.out.tsv
    hits = SHIGATYPER_MODULE.out.hits
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
