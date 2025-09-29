//
// btyper3 - Taxonomic classification of Bacillus cereus group isolates
//
include { BTYPER3 as BTYPER3_MODULE } from '../../modules/btyper3/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BTYPER3 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    BTYPER3_MODULE(fasta)

    // Merge results
    BTYPER3_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'btyper3'], tsv]}.set{ ch_merge_btyper3 }
    CSVTK_CONCAT(ch_merge_btyper3, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = BTYPER3_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = BTYPER3_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = BTYPER3_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = BTYPER3_MODULE.out.nf_begin.mix(
        BTYPER3_MODULE.out.nf_err,
        BTYPER3_MODULE.out.nf_log,
        BTYPER3_MODULE.out.nf_out,
        BTYPER3_MODULE.out.nf_run,
        BTYPER3_MODULE.out.nf_sh,
        BTYPER3_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = BTYPER3_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
