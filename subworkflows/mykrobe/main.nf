//
// mykrobe - Antimicrobial resistance detection for specific species
//
include { MYKROBE_PREDICT }  from '../../modules/mykrobe/predict/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MYKROBE {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    MYKROBE_PREDICT(reads, params.mykrobe_species)
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())
    ch_logs = ch_logs.mix(MYKROBE_PREDICT.out.logs)
    
    // Merge results
    MYKROBE_PREDICT.out.csv.collect{_meta, csv -> csv}.map{ csv -> [[id:'mykrobe'], csv]}.set{ ch_merge_mykrobe }
    CSVTK_CONCAT(ch_merge_mykrobe, 'csv', 'csv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    csv = MYKROBE_PREDICT.out.csv
    json = MYKROBE_PREDICT.out.json
    merged_csv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = MYKROBE_PREDICT.out.nf_begin.mix(
        MYKROBE_PREDICT.out.nf_err,
        MYKROBE_PREDICT.out.nf_log,
        MYKROBE_PREDICT.out.nf_out,
        MYKROBE_PREDICT.out.nf_run,
        MYKROBE_PREDICT.out.nf_sh,
        MYKROBE_PREDICT.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions
}
