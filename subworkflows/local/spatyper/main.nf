//
// SpaTyper - Computational method for finding spa types in Staphylococcus aureus
//
REPEATS = params.repeats ? file(params.repeats, checkIfExists: true) : []
REPEAT_ORDER = params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []
spatyper_args = params.do_enrich ? '--do_enrich' : ''

include { SPATYPER as SPATYPER_MODULE } from '../../../modules/nf-core/modules/spatyper/main' addParams( options: [args: "${spatyper_args}", is_module: true] )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow SPATYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    SPATYPER_MODULE(fasta, REPEATS, REPEAT_ORDER)
    ch_versions = ch_versions.mix(SPATYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        SPATYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'spatyper'], tsv]}.set{ ch_merge_spatyper }
        CSVTK_CONCAT(ch_merge_spatyper, 'tsv', 'tsv')
    }

    emit:
    spatyper_tsv = SPATYPER_MODULE.out.tsv
    spatyper_all_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
