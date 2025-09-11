//
// mlst - Automatic MLST calling from assembled contigs
//
include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MLST {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ mlst_db ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    if (params.mlst_db) {
        MLST_MODULE(fasta, file(params.mlst_db))
    } else {
        MLST_MODULE(fasta, db)
    }

    ch_versions = ch_versions.mix(MLST_MODULE.out.versions)
    ch_logs = ch_logs.mix(MLST_MODULE.out.logs)

    // Merge results
    MLST_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'mlst'], tsv]}.set{ ch_merge_mlst }
    CSVTK_CONCAT(ch_merge_mlst, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = MLST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        MLST_MODULE.out.nf_begin,
        MLST_MODULE.out.nf_err,
        MLST_MODULE.out.nf_log,
        MLST_MODULE.out.nf_out,
        MLST_MODULE.out.nf_run,
        MLST_MODULE.out.nf_sh,
        MLST_MODULE.out.nf_trace
    )
    versions = ch_versions
}
