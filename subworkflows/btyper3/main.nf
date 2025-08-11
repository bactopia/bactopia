//
// btyper3 - Taxonomic classification of Bacillus cereus group isolates
//
include { BTYPER3 as BTYPER3_MODULE } from '../../modules/btyper3/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BTYPER3 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    BTYPER3_MODULE(fasta)
    ch_versions = ch_versions.mix(BTYPER3_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(BTYPER3_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BTYPER3_MODULE.out.nf_logs)

    // Merge results
    BTYPER3_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'btyper3'], tsv] }.set{ ch_merge_btyper3 }
    CSVTK_CONCAT(ch_merge_btyper3, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = BTYPER3_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
