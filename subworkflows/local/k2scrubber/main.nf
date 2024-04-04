//
// k2scrubber - Scrub human reads from FASTQ files using Kraken2 and human pangenome reference
//

include { initOptions } from '../../../lib/nf/functions'

// HPRC download options
wget_opts = initOptions(params.containsKey("options") ? params.options : [:], 'wget')
wget_opts.args = [
    params.no_check_certificate ? "--no-check-certificate" : "",
    "-O k2_HPRC_20230810.tar.gz",
    "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1",
].join(' ').replaceAll("\\s{2,}", " ").trim()
wget_opts.args2 = "k2_HPRC_20230810.tar.gz"

// Kraken2 Options
k2_options = initOptions(params.containsKey("options") ? params.options : [:], 'kraken2')
k2_options.args = [
    params.quick_mode ? "--quick" : "",
    params.use_mpa_style ? "--use-mpa-style" : "",
    params.report_zero_counts ? "--report-zero-counts" : "",
    params.report_minimizer_data ? "--report-minimizer-data" : "",
    params.use_names ? "--use-names" : "",
    params.memory_mapping ? "--memory-mapping" : "",
    "--confidence ${params.confidence}",
    "--minimum-hit-groups ${params.minimum_hit_groups}",
    "--minimum-base-quality ${params.minimum_base_quality}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { CUSTOM_WGET as WGET_HPRC } from '../../../modules/nf-core/custom/wget/main' addParams( options: wget_opts )
include { KRAKEN2 } from '../../../modules/nf-core/kraken2/main' addParams( options: k2_options )

workflow K2SCRUBBER {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    WGET_HPRC()
    KRAKEN2(reads, WGET_HPRC.out.download)
    ch_versions = ch_versions.mix(KRAKEN2.out.versions)

    emit:
    scrubbed = KRAKEN2.out.unclassified
    scrubbed_extra = KRAKEN2.out.unclassified_extra
    human = KRAKEN2.out.classified
    kraken2_report = KRAKEN2.out.kraken2_report
    scrub_report = KRAKEN2.out.scrub_report
    versions = ch_versions // channel: [ versions.yml ]
}
