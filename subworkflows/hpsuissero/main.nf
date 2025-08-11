//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow HPSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    HPSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(HPSUISSERO_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(HPSUISSERO_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(HPSUISSERO_MODULE.out.nf_logs)

    // Merge results
    HPSUISSERO_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'hpsuissero'], tsv] }.set{ ch_merge_hpsuissero }
    CSVTK_CONCAT(ch_merge_hpsuissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = HPSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
