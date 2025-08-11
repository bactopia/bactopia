//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    LISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(LISSERO_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(LISSERO_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(LISSERO_MODULE.out.nf_logs)

    // Merge results
    LISSERO_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'lissero'], tsv] }.set{ ch_merge_lissero }
    CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
