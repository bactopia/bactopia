//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'hpsuissero')
options.is_module = params.wf == 'hpsuissero' ? true : false

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../../modules/nf-core/hpsuissero/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'hpsuisero'] )
}

workflow HPSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_hpsuissero = Channel.empty()

    HPSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(HPSUISSERO_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        HPSUISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'hpsuissero'], tsv]}.set{ ch_merge_hpsuissero }
        CSVTK_CONCAT(ch_merge_hpsuissero, 'tsv', 'tsv')
        ch_merged_hpsuissero = ch_merged_hpsuissero.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = HPSUISSERO_MODULE.out.tsv
    merged_tsv = ch_merged_hpsuissero
    versions = ch_versions
}
