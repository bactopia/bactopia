//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'lissero')
options.is_module = params.wf == 'lissero' ? true : false
options.args = [
    "--min_id ${params.min_id}",
    "--min_cov ${params.min_cov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { LISSERO as LISSERO_MODULE } from '../../../modules/nf-core/modules/lissero/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'lissero'] )
}

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_lissero = Channel.empty()

    LISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(LISSERO_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        LISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'lissero'], tsv]}.set{ ch_merge_lissero }
        CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')
        ch_merged_lissero = ch_merged_lissero.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = ch_merged_lissero
    versions = ch_versions
}
