//
// lissero - Serogroup typing prediction for Listeria monocytogenes
//
include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow LISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    LISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(LISSERO_MODULE.out.versions.first()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions))
    ch_logs = ch_logs.mix(LISSERO_MODULE.out.logs)

    // Merge results
    LISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'lissero'], tsv]}.set{ ch_merge_lissero }
    CSVTK_CONCAT(ch_merge_lissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = LISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
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
    versions = ch_versions
}
