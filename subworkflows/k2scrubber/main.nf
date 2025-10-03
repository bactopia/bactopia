//
// k2scrubber - Scrub human reads from FASTQ files using Kraken2 and human pangenome reference
//
include { WGET    } from '../../modules/wget/main'
include { KRAKEN2 } from '../../modules/kraken2/main'

workflow K2SCRUBBER {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    WGET([
        "name": "k2scrubber",
        "save_as": "k2_HPRC_20230810.tar.gz",
        "url": "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1"
    ])
    KRAKEN2(reads, WGET.out.download)

    emit:
    // Individual outputs
    human = KRAKEN2.out.classified
    kraken2_report = KRAKEN2.out.kraken2_report
    scrub_report = KRAKEN2.out.scrub_report
    scrub_special_report = KRAKEN2.out.scrub_special_report
    scrubbed = KRAKEN2.out.unclassified
    scrubbed_extra = KRAKEN2.out.unclassified_extra

    // Generic aggregate outputs
    results = KRAKEN2.out.classified.mix(
        KRAKEN2.out.kraken2_report,
        KRAKEN2.out.scrub_report,
        KRAKEN2.out.unclassified,
    )
    logs = KRAKEN2.out.logs
    nf_logs = KRAKEN2.out.nf_begin.mix(
        KRAKEN2.out.nf_err,
        KRAKEN2.out.nf_log,
        KRAKEN2.out.nf_out,
        KRAKEN2.out.nf_run,
        KRAKEN2.out.nf_sh,
        KRAKEN2.out.nf_trace
    )
    versions = KRAKEN2.out.versions
}
