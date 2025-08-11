//
// kraken2 - Taxonomic classification of sequence reads 
//
include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'

workflow KRAKEN2 {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    DATABASE = params.kraken2_db ? file(params.kraken2_db) : []

    KRAKEN2_MODULE(reads, DATABASE)
    ch_versions = ch_versions.mix(KRAKEN2_MODULE.out.versions)
    ch_logs = ch_logs.mix(KRAKEN2_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KRAKEN2_MODULE.out.nf_logs)

    emit:
    classified = KRAKEN2_MODULE.out.classified
    unclassified = KRAKEN2_MODULE.out.unclassified
    kraken2_report = KRAKEN2_MODULE.out.kraken2_report
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
