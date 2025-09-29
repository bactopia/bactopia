//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    LISSERO_MODULE(fasta)

    // Merge results
    LISSERO_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'lissero'], tsv]}.set{ ch_merge_lissero }
    CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = LISSERO_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = LISSERO_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        LISSERO_MODULE.out.nf_begin,
        LISSERO_MODULE.out.nf_err,
        LISSERO_MODULE.out.nf_log,
        LISSERO_MODULE.out.nf_out,
        LISSERO_MODULE.out.nf_run,
        LISSERO_MODULE.out.nf_sh,
        LISSERO_MODULE.out.nf_trace
    )
    versions = LISSERO_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
