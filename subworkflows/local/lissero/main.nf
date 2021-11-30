//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
lissero_opts = [
    "--min_id ${params.min_id}",
    "--min_cov ${params.min_cov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { LISSERO as LISSERO_MODULE } from '../../../modules/nf-core/modules/lissero/main' addParams( options: [ args: "${lissero_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    LISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(LISSERO_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        LISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'lissero'], tsv]}.set{ ch_merge_lissero }
        CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }

    emit:
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
