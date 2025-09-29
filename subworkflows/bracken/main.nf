//
// bracken - Estimate taxonomic abundance of samples from Kraken2 results
//
include { BRACKEN as BRACKEN_MODULE } from '../../modules/bracken/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../modules/csvtk/concat/main'

workflow BRACKEN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    database

    main:
    BRACKEN_MODULE(reads, database)

    // Merge Bracken Primary/Secondary Species abundance
    BRACKEN_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'bracken-species-abundance'], tsv]}.set{ ch_merge_tsv }
    CSVTK_CONCAT_TSV(ch_merge_tsv, 'tsv', 'tsv')

    // Merge Bracken adjusted abundance
    BRACKEN_MODULE.out.adjusted_abundances.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'bracken-adjusted'], tsv]}.set{ ch_merge_adjusted_abundances }
    CSVTK_CONCAT_ADJUSTED(ch_merge_adjusted_abundances, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BRACKEN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT_TSV.out.csv
    abundances = BRACKEN_MODULE.out.abundances
    adjusted_abundances = BRACKEN_MODULE.out.adjusted_abundances
    bracken_report = BRACKEN_MODULE.out.bracken_report
    classification = BRACKEN_MODULE.out.classification
    classified = BRACKEN_MODULE.out.classified
    kraken2_output = BRACKEN_MODULE.out.kraken2_output
    kraken2_report = BRACKEN_MODULE.out.kraken2_report
    merged_adjusted_abundances = CSVTK_CONCAT_ADJUSTED.out.csv
    unclassified = BRACKEN_MODULE.out.unclassified

    // Generic aggregate outputs
    results = BRACKEN_MODULE.out.tsv.mix(
        CSVTK_CONCAT_TSV.out.csv,
        BRACKEN_MODULE.out.abundances,
        BRACKEN_MODULE.out.adjusted_abundances,
        BRACKEN_MODULE.out.bracken_report,
        BRACKEN_MODULE.out.classification,
        BRACKEN_MODULE.out.classified,
        BRACKEN_MODULE.out.kraken2_output,
        BRACKEN_MODULE.out.kraken2_report,
        CSVTK_CONCAT_ADJUSTED.out.csv,
        BRACKEN_MODULE.out.unclassified
    )
    logs = BRACKEN_MODULE.out.logs.mix(
        CSVTK_CONCAT_TSV.out.logs,
        CSVTK_CONCAT_ADJUSTED.out.logs
    )
    nf_logs = BRACKEN_MODULE.out.nf_begin.mix(
        BRACKEN_MODULE.out.nf_err,
        BRACKEN_MODULE.out.nf_log,
        BRACKEN_MODULE.out.nf_out,
        BRACKEN_MODULE.out.nf_run,
        BRACKEN_MODULE.out.nf_sh,
        BRACKEN_MODULE.out.nf_trace
    )
    versions = BRACKEN_MODULE.out.versions.mix(
        CSVTK_CONCAT_TSV.out.versions,
        CSVTK_CONCAT_ADJUSTED.out.versions
    )
}
