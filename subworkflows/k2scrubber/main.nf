//
// k2scrubber - Scrub human reads from FASTQ files using Kraken2 and human pangenome reference
//
include { CUSTOM_WGET as WGET_HPRC } from '../../modules/custom/wget/main'
include { KRAKEN2 } from '../../modules/kraken2/main'

workflow K2SCRUBBER {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    WGET_HPRC()
    KRAKEN2(reads, WGET_HPRC.out.download)

    emit:
    // Individual outputs
    human = KRAKEN2.out.classified
    kraken2_report = KRAKEN2.out.kraken2_report
    scrub_report = KRAKEN2.out.scrub_report
    scrubbed = KRAKEN2.out.unclassified
    scrubbed_extra = KRAKEN2.out.unclassified_extra

    // Generic aggregate outputs
    results = KRAKEN2.out.classified.mix(
        KRAKEN2.out.kraken2_report,
        KRAKEN2.out.scrub_report,
        KRAKEN2.out.unclassified,
    )
    logs = KRAKEN2.out.logs.mix(
        WGET_HPRC.out.logs
    )
    nf_logs = KRAKEN2.out.nf_begin.mix(
        KRAKEN2.out.nf_err,
        KRAKEN2.out.nf_log,
        KRAKEN2.out.nf_out,
        KRAKEN2.out.nf_run,
        KRAKEN2.out.nf_sh,
        KRAKEN2.out.nf_trace,
        WGET_HPRC.out.nf_begin,
        WGET_HPRC.out.nf_err,
        WGET_HPRC.out.nf_log,
        WGET_HPRC.out.nf_out,
        WGET_HPRC.out.nf_run,
        WGET_HPRC.out.nf_sh,
        WGET_HPRC.out.nf_trace
    )
    versions = WGET_HPRC.out.versions.mix(
        KRAKEN2.out.versions
    )
}
