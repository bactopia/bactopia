//
// mlst - Automatic MLST calling from assembled contigs
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mlst')
options.is_module = params.wf == 'mlst' ? true : false
options.args = [
    params.nopath ? "--nopath" : "",
    params.scheme ? "--scheme ${params.scheme}" : "",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}",
    "--minscore ${params.minscore}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MLST as MLST_MODULE } from '../../../modules/nf-core/modules/mlst/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [args: '--no-header-row', publish_to_base: true, logs_subdir: options.is_module ? '' : 'mlst'] )

}
workflow MLST {
    take:
    fasta // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_mlst = Channel.empty()

    MLST_MODULE(fasta)
    ch_versions = ch_versions.mix(MLST_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        MLST_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'mlst'], tsv]}.set{ ch_merge_mlst }
        CSVTK_CONCAT(ch_merge_mlst, 'tsv', 'tsv')
        ch_merged_mlst = ch_merged_mlst.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = MLST_MODULE.out.tsv
    merged_tsv = ch_merged_mlst
    versions = ch_versions // channel: [ versions.yml ]
}
