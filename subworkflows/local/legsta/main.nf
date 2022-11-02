//
// legsta - Typing of Legionella pneumophila assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'legsta')
options.is_module = params.wf == 'legsta' ? true : false
options.args = params.noheader ? '--noheader' : ''

include { LEGSTA as LEGSTA_MODULE } from '../../../modules/nf-core/legsta/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'legsta'] )
}

workflow LEGSTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_legsta = Channel.empty()

    LEGSTA_MODULE(fasta)
    ch_versions = ch_versions.mix(LEGSTA_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        LEGSTA_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'legsta'], tsv]}.set{ ch_merge_legsta }
        CSVTK_CONCAT(ch_merge_legsta, 'tsv', 'tsv')
        ch_merged_legsta = ch_merged_legsta.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = LEGSTA_MODULE.out.tsv
    merged_tsv = ch_merged_legsta
    versions = ch_versions
}
