//
// midas - Estimate species abundances from FASTQ files
//
include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MIDAS {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    DATABASE = params.midas_db ? file(params.midas_db) : []

    MIDAS_SPECIES(reads, DATABASE)
    ch_versions = ch_versions.mix(MIDAS_SPECIES.out.versions)
    ch_logs = ch_logs.mix(MIDAS_SPECIES.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MIDAS_SPECIES.out.nf_logs)

    // Merge results
    MIDAS_SPECIES.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'midas'], tsv] }.set{ ch_merge_midas }
    CSVTK_CONCAT(ch_merge_midas, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = MIDAS_SPECIES.out.tsv
    abundances = MIDAS_SPECIES.out.abundances
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
