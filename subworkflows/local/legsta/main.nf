//
// legsta - Typing of Legionella pneumophila assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'legsta')
options.args = params.noheader ? '--noheader' : ''

include { LEGSTA as LEGSTA_MODULE } from '../../../modules/nf-core/legsta/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'legsta'] )

workflow LEGSTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    LEGSTA_MODULE(fasta)
    ch_versions = ch_versions.mix(LEGSTA_MODULE.out.versions.first())

    // Merge results
    LEGSTA_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'legsta'], tsv]}.set{ ch_merge_legsta }
    CSVTK_CONCAT(ch_merge_legsta, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = LEGSTA_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
