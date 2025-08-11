//
// clermontyping - in silico phylotyping of Escherichia genus
//
include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../modules/clermontyping/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CLERMONTYPING {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    CLERMONTYPING_MODULE(fasta)
    ch_versions = ch_versions.mix(CLERMONTYPING_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(CLERMONTYPING_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CLERMONTYPING_MODULE.out.nf_logs)

    // Merge results
    CLERMONTYPING_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'clermontyping'], tsv] }.set{ ch_merge_clermontyping }
    CSVTK_CONCAT(ch_merge_clermontyping, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = CLERMONTYPING_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
