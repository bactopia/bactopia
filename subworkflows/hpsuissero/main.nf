//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow HPSUISSERO {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    HPSUISSERO_MODULE(fasta)
    ch_versions = ch_versions.mix(HPSUISSERO_MODULE.out.versions.first()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions))
    ch_logs = ch_logs.mix(HPSUISSERO_MODULE.out.logs)

    // Merge results
    HPSUISSERO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'hpsuissero'], tsv]}.set{ ch_merge_hpsuissero }
    CSVTK_CONCAT(ch_merge_hpsuissero, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = HPSUISSERO_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        HPSUISSERO_MODULE.out.nf_begin,
        HPSUISSERO_MODULE.out.nf_err,
        HPSUISSERO_MODULE.out.nf_log,
        HPSUISSERO_MODULE.out.nf_out,
        HPSUISSERO_MODULE.out.nf_run,
        HPSUISSERO_MODULE.out.nf_sh,
        HPSUISSERO_MODULE.out.nf_trace
    )
    versions = ch_versions
}
