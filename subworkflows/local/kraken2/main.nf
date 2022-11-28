//
// kraken2 - Taxonomic classification of sequence reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'kraken2')
options.is_module = params.wf == 'kraken2' ? true : false
options.args = [
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
DATABASE = params.kraken2_db ? file(params.kraken2_db) : []

include { KRAKEN2 as KRAKEN2_MODULE } from '../../../modules/nf-core/kraken2/main' addParams( options: options )

workflow KRAKEN2 {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    KRAKEN2_MODULE(reads, DATABASE)
    ch_versions = ch_versions.mix(KRAKEN2_MODULE.out.versions)

    emit:
    classified = KRAKEN2_MODULE.out.classified
    unclassified = KRAKEN2_MODULE.out.unclassified
    report = KRAKEN2_MODULE.out.report
    versions = ch_versions // channel: [ versions.yml ]
}
