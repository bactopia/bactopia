//
// legsta - Typing of Legionella pneumophila assemblies
//
include { LEGSTA as LEGSTA_MODULE } from '../../modules/legsta/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LEGSTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    LEGSTA_MODULE(fasta)
    ch_versions = ch_versions.mix(LEGSTA_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(LEGSTA_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(LEGSTA_MODULE.out.nf_logs)

    // Merge results
    LEGSTA_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'legsta'], tsv] }.set{ ch_merge_legsta }
    CSVTK_CONCAT(ch_merge_legsta, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = LEGSTA_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
