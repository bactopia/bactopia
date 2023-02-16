//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'hpsuissero')

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../../modules/nf-core/hpsuissero/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'hpsuisero'] )

workflow HPSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    HPSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(HPSUISSERO_MODULE.out.versions.first())

    // Merge results
    HPSUISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'hpsuissero'], tsv]}.set{ ch_merge_hpsuissero }
    CSVTK_CONCAT(ch_merge_hpsuissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = HPSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
