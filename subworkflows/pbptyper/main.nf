//
// pbptyper - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae
//
include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PBPTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    PBPTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(PBPTYPER_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(PBPTYPER_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PBPTYPER_MODULE.out.nf_logs)

    // Merge results
    PBPTYPER_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'pbptyper'], tsv] }.set{ ch_merge_pbptyper }
    CSVTK_CONCAT(ch_merge_pbptyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = PBPTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PBPTYPER_MODULE.out.blast
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
