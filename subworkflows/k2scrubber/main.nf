//
// k2scrubber - Scrub human reads from FASTQ files using Kraken2 and human pangenome reference
//
nextflow.preview.types = true

include { WGET         } from '../../modules/wget/main'
include { KRAKEN2      } from '../../modules/kraken2/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow K2SCRUBBER {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    WGET([
        "name": "k2scrubber",
        "save_as": "k2_HPRC_20230810.tar.gz",
        "url": "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1"
    ])
    KRAKEN2(reads, WGET.out.download)

    emit:
    // Individual outputs
    human: Channel<Tuple<Map, Path>> = KRAKEN2.out.classified
    kraken2_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.kraken2_report
    scrub_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_special_report
    scrubbed: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified
    scrubbed_extra: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified_extra

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2.out.classified,
        KRAKEN2.out.kraken2_report,
        KRAKEN2.out.scrub_report,
        KRAKEN2.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = KRAKEN2.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = KRAKEN2.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = KRAKEN2.out.versions
}
