//
// ssuissero - Serotype prediction of Streptococcus suis assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ssuissero')
options.is_module = params.wf == 'ssuissero' ? true : false

include { SSUISSERO as SSUISSERO_MODULE } from '../../../modules/nf-core/modules/ssuissero/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'ssuissero'] )
}

workflow SSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_ssuissero = Channel.empty()

    SSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(SSUISSERO_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SSUISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ssuissero'], tsv]}.set{ ch_merge_ssuissero }
        CSVTK_CONCAT(ch_merge_ssuissero, 'tsv', 'tsv')
        ch_merged_ssuissero = ch_merged_ssuissero.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = SSUISSERO_MODULE.out.tsv
    merged_tsv = ch_merged_ssuissero
    versions = ch_versions
}
