//
// SpaTyper - Computational method for finding spa types in Staphylococcus aureus
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'spatyper')
options.args = params.do_enrich ? '--do_enrich' : ''
REPEATS = params.repeats ? file(params.repeats, checkIfExists: true) : []
REPEAT_ORDER = params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []

include { SPATYPER as SPATYPER_MODULE } from '../../../modules/nf-core/spatyper/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'spatyper-concat', process_name: params.merge_folder] )

workflow SPATYPER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    SPATYPER_MODULE(fasta, REPEATS, REPEAT_ORDER)
    ch_versions = ch_versions.mix(SPATYPER_MODULE.out.versions.first())

    // Merge results
    SPATYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'spatyper'], tsv]}.set{ ch_merge_spatyper }
    CSVTK_CONCAT(ch_merge_spatyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SPATYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
