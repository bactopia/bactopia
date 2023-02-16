//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'lissero')
options.args = [
    "--min_id ${params.min_id}",
    "--min_cov ${params.min_cov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { LISSERO as LISSERO_MODULE } from '../../../modules/nf-core/lissero/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'lissero'] )

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    LISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(LISSERO_MODULE.out.versions.first())

    // Merge results
    LISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'lissero'], tsv]}.set{ ch_merge_lissero }
    CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
