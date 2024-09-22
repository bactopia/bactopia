//
// bracken - Estimate taxonomic abundance of samples from Kraken2 results
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'bracken')
options.is_module = params.wf == 'bracken' ? true : false
options.args = [
    params.kraken2_quick_mode ? "--quick" : "",
    params.kraken2_use_mpa_style ? "--use-mpa-style" : "",
    params.kraken2_report_zero_counts ? "--report-zero-counts" : "",
    params.kraken2_report_minimizer_data ? "--report-minimizer-data" : "",
    params.kraken2_use_names ? "--use-names" : "",
    params.kraken2_memory_mapping ? "--memory-mapping" : "",
    "--confidence ${params.kraken2_confidence}",
    "--minimum-hit-groups ${params.kraken2_minimum_hit_groups}",
    "--minimum-base-quality ${params.kraken2_minimum_base_quality}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

options.args2 = [
    "-l ${params.bracken_level}",
    "-t ${params.bracken_threshold}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE = params.kraken2_db ? file(params.kraken2_db) : []

include { BRACKEN as BRACKEN_MODULE } from '../../../modules/nf-core/bracken/main' addParams( options: options )
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'bracken-species-abundance-concat', process_name: params.merge_folder])
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'bracken-adjusted-concat', process_name: params.merge_folder])

workflow BRACKEN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    BRACKEN_MODULE(reads, DATABASE)
    ch_versions = ch_versions.mix(BRACKEN_MODULE.out.versions)

    // Merge Bracken Primary/Secondary Species abundance
    BRACKEN_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'bracken-species-abundance'], tsv]}.set{ ch_merge_tsv }
    CSVTK_CONCAT_TSV(ch_merge_tsv, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_TSV.out.versions)

    // Merge Bracken adjusted abundance
    BRACKEN_MODULE.out.adjusted_abundances.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'bracken-adjusted'], tsv]}.set{ ch_merge_adjusted_abundances }
    CSVTK_CONCAT_ADJUSTED(ch_merge_adjusted_abundances, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_ADJUSTED.out.versions)

    emit:
    tsv = BRACKEN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT_TSV.out.csv
    classified = BRACKEN_MODULE.out.classified
    unclassified = BRACKEN_MODULE.out.unclassified
    kraken2_report = BRACKEN_MODULE.out.kraken2_report
    kraken2_output = BRACKEN_MODULE.out.kraken2_output
    bracken_report = BRACKEN_MODULE.out.bracken_report
    abundances = BRACKEN_MODULE.out.abundances
    adjusted_abundances = BRACKEN_MODULE.out.adjusted_abundances
    merged_adjusted_abundances = CSVTK_CONCAT_ADJUSTED.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
