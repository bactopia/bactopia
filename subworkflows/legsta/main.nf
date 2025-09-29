//
// legsta - Typing of Legionella pneumophila assemblies
//
include { LEGSTA as LEGSTA_MODULE } from '../../modules/legsta/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LEGSTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    LEGSTA_MODULE(fasta)

    // Merge results
    LEGSTA_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'legsta'], tsv]}.set{ ch_merge_legsta }
    CSVTK_CONCAT(ch_merge_legsta, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = LEGSTA_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = LEGSTA_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = LEGSTA_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        LEGSTA_MODULE.out.nf_begin,
        LEGSTA_MODULE.out.nf_err,
        LEGSTA_MODULE.out.nf_log,
        LEGSTA_MODULE.out.nf_out,
        LEGSTA_MODULE.out.nf_run,
        LEGSTA_MODULE.out.nf_sh,
        LEGSTA_MODULE.out.nf_trace
    )
    versions = LEGSTA_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
