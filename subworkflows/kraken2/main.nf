//
// kraken2 - Taxonomic classification of sequence reads 
//
include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'

workflow KRAKEN2 {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    database

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    KRAKEN2_MODULE(reads, database)
    ch_versions = ch_versions.mix(KRAKEN2_MODULE.out.versions)
    ch_logs = ch_logs.mix(KRAKEN2_MODULE.out.logs)

    emit:
    classified = KRAKEN2_MODULE.out.classified
    kraken2_report = KRAKEN2_MODULE.out.kraken2_report
    unclassified = KRAKEN2_MODULE.out.unclassified
    logs = ch_logs
    nf_logs = KRAKEN2_MODULE.out.nf_begin.mix(
        KRAKEN2_MODULE.out.nf_err,
        KRAKEN2_MODULE.out.nf_log,
        KRAKEN2_MODULE.out.nf_out,
        KRAKEN2_MODULE.out.nf_run,
        KRAKEN2_MODULE.out.nf_sh,
        KRAKEN2_MODULE.out.nf_trace
    )
    versions = ch_versions
}
