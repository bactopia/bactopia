//
// meningotype - Serotyping of Neisseria meningitidis
//
include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    MENINGOTYPE_MODULE(fasta)
    ch_versions = ch_versions.mix(MENINGOTYPE_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(MENINGOTYPE_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MENINGOTYPE_MODULE.out.nf_logs)

    // Merge results
    MENINGOTYPE_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'meningotype'], tsv] }.set{ ch_merge_meningotype }
    CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
