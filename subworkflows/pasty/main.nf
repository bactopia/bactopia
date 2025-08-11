//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    PASTY_MODULE(fasta)
    ch_versions = ch_versions.mix(PASTY_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(PASTY_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PASTY_MODULE.out.nf_logs)

    // Merge results
    PASTY_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'pasty'], tsv] }.set{ ch_merge_pasty }
    CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PASTY_MODULE.out.blast
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
