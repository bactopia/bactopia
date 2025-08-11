//
// ectyper - In-silico prediction of Escherichia coli serotype
//
include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    ECTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(ECTYPER_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(ECTYPER_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ECTYPER_MODULE.out.nf_logs)

    // Merge results
    ECTYPER_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'ectyper'], tsv] }.set{ ch_merge_ectyper }
    CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = ECTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    txt = ECTYPER_MODULE.out.txt
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
