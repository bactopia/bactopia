//
// kraken2 - Taxonomic classification of sequence reads 
//
nextflow.preview.types = true

include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'

workflow KRAKEN2 {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    database

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    // Individual outputs
    classified = KRAKEN2_MODULE.out.classified
    kraken2_report = KRAKEN2_MODULE.out.kraken2_report
    unclassified = KRAKEN2_MODULE.out.unclassified

    // Generic aggregate outputs
    results = KRAKEN2_MODULE.out.classified.mix(
        KRAKEN2_MODULE.out.kraken2_report,
        KRAKEN2_MODULE.out.unclassified
    )
    logs = KRAKEN2_MODULE.out.logs
    nf_logs = KRAKEN2_MODULE.out.nf_begin.mix(
        KRAKEN2_MODULE.out.nf_err,
        KRAKEN2_MODULE.out.nf_log,
        KRAKEN2_MODULE.out.nf_out,
        KRAKEN2_MODULE.out.nf_run,
        KRAKEN2_MODULE.out.nf_sh,
        KRAKEN2_MODULE.out.nf_trace
    )
    versions = KRAKEN2_MODULE.out.versions
}
