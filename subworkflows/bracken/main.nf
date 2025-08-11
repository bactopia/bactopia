//
// bracken - Estimate taxonomic abundance of samples from Kraken2 results
//
include { BRACKEN as BRACKEN_MODULE } from '../../modules/bracken/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../modules/csvtk/concat/main'

workflow BRACKEN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    DATABASE = params.kraken2_db ? file(params.kraken2_db) : []

    BRACKEN_MODULE(reads, DATABASE)
    ch_versions = ch_versions.mix(BRACKEN_MODULE.out.versions)
    ch_logs = ch_logs.mix(BRACKEN_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BRACKEN_MODULE.out.nf_logs)

    // Merge Bracken Primary/Secondary Species abundance
    BRACKEN_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'bracken-species-abundance'], tsv] }.set{ ch_merge_tsv }
    CSVTK_CONCAT_TSV(ch_merge_tsv, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_TSV.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT_TSV.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT_TSV.out.nf_logs)

    // Merge Bracken adjusted abundance
    BRACKEN_MODULE.out.adjusted_abundances.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'bracken-adjusted'], tsv] }.set{ ch_merge_adjusted_abundances }
    CSVTK_CONCAT_ADJUSTED(ch_merge_adjusted_abundances, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT_ADJUSTED.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT_ADJUSTED.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT_ADJUSTED.out.nf_logs)

    emit:
    tsv = BRACKEN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT_TSV.out.csv
    classified = BRACKEN_MODULE.out.classified
    unclassified = BRACKEN_MODULE.out.unclassified
    kraken2_report = BRACKEN_MODULE.out.kraken2_report
    kraken2_output = BRACKEN_MODULE.out.kraken2_output
    bracken_report = BRACKEN_MODULE.out.bracken_report
    abundances = BRACKEN_MODULE.out.abundances
    classification = BRACKEN_MODULE.out.classification
    adjusted_abundances = BRACKEN_MODULE.out.adjusted_abundances
    merged_adjusted_abundances = CSVTK_CONCAT_ADJUSTED.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
