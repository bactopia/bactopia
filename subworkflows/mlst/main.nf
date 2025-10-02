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
    MLST_MODULE(fasta, db)

    // Merge results
    MLST_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'mlst'], tsv]}.set{ ch_merge_mlst }
    CSVTK_CONCAT(ch_merge_mlst, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = MLST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = MLST_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = MLST_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = MLST_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
