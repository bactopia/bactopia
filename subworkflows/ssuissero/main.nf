//
// ssuissero - Serotype prediction of Streptococcus suis assemblies
//
include { SSUISSERO as SSUISSERO_MODULE } from '../../modules/ssuissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    SSUISSERO_MODULE(fasta)

    // Merge results
    SSUISSERO_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'ssuissero'], tsv]}.set{ ch_merge_ssuissero }
    CSVTK_CONCAT(ch_merge_ssuissero, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = SSUISSERO_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = SSUISSERO_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SSUISSERO_MODULE.out.nf_begin.mix(
        SSUISSERO_MODULE.out.nf_err,
        SSUISSERO_MODULE.out.nf_log,
        SSUISSERO_MODULE.out.nf_out,
        SSUISSERO_MODULE.out.nf_run,
        SSUISSERO_MODULE.out.nf_sh,
        SSUISSERO_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SSUISSERO_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
