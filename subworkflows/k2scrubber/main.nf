//
// k2scrubber - Scrub human reads from FASTQ files using Kraken2 and human pangenome reference
//
include { CUSTOM_WGET as WGET_HPRC } from '../../modules/custom/wget/main'
include { KRAKEN2 } from '../../modules/kraken2/main'

workflow K2SCRUBBER {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    WGET_HPRC()
    ch_versions = ch_versions.mix(WGET_HPRC.out.versions)
    ch_logs = ch_logs.mix(WGET_HPRC.out.logs)
    ch_nf_logs = ch_nf_logs.mix(WGET_HPRC.out.nf_logs)

    KRAKEN2(reads, WGET_HPRC.out.download)
    ch_versions = ch_versions.mix(KRAKEN2.out.versions)
    ch_logs = ch_logs.mix(KRAKEN2.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KRAKEN2.out.nf_logs)

    emit:
    scrubbed = KRAKEN2.out.unclassified
    scrubbed_extra = KRAKEN2.out.unclassified_extra
    human = KRAKEN2.out.classified
    kraken2_report = KRAKEN2.out.kraken2_report
    scrub_report = KRAKEN2.out.scrub_report
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
