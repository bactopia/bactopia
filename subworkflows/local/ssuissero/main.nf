//
// ssuissero - Serotype prediction of Streptococcus suis assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ssuissero')

include { SSUISSERO as SSUISSERO_MODULE } from '../../../modules/nf-core/ssuissero/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'ssuissero'] )

workflow SSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(SSUISSERO_MODULE.out.versions.first())

    // Merge results
    SSUISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ssuissero'], tsv]}.set{ ch_merge_ssuissero }
    CSVTK_CONCAT(ch_merge_ssuissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
