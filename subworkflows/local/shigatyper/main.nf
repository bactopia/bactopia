//
// shigatyper - Shigella serotype from Illumina or Oxford Nanopore reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'shigatyper')
options.is_module = params.wf == 'shigatyper' ? true : false
include { SHIGATYPER as SHIGATYPER_MODULE } from '../../../modules/nf-core/shigatyper/main' addParams( options: options )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'shigatyper'] )
}

workflow SHIGATYPER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_shigatyper = Channel.empty()

    SHIGATYPER_MODULE(reads)
    ch_versions = ch_versions.mix(SHIGATYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SHIGATYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'shigatyper'], tsv]}.set{ ch_merge_shigatyper }
        CSVTK_CONCAT(ch_merge_shigatyper, 'tsv', 'tsv')
        ch_merged_shigatyper = ch_merged_shigatyper.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SHIGATYPER_MODULE.out.tsv
    hits = SHIGATYPER_MODULE.out.hits
    merged_tsv = ch_merged_shigatyper
    versions = ch_versions
}
