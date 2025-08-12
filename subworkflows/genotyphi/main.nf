//
// genotyphi - Salmonella Typhi genotyping with Mykrobe outputs
//
include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { GENOTYPHI_PARSE  } from '../../modules/genotyphi/parse/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow GENOTYPHI {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    MYKROBE_PREDICT(reads, "typhi")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())
    ch_logs = ch_logs.mix(MYKROBE_PREDICT.out.logs)
    
    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    ch_versions = ch_versions.mix(GENOTYPHI_PARSE.out.versions.first())
    ch_logs = ch_logs.mix(GENOTYPHI_PARSE.out.logs)
    
    // Merge results
    GENOTYPHI_PARSE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'genotyphi'], tsv]}.set{ ch_merge_genotyphi }
    CSVTK_CONCAT(ch_merge_genotyphi, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = GENOTYPHI_PARSE.out.tsv
    csv = MYKROBE_PREDICT.out.csv
    json = MYKROBE_PREDICT.out.json
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = MYKROBE_PREDICT.out.nf_begin.mix(
        MYKROBE_PREDICT.out.nf_err,
        MYKROBE_PREDICT.out.nf_log,
        MYKROBE_PREDICT.out.nf_out,
        MYKROBE_PREDICT.out.nf_run,
        MYKROBE_PREDICT.out.nf_sh,
        MYKROBE_PREDICT.out.nf_trace,
        GENOTYPHI_PARSE.out.nf_begin,
        GENOTYPHI_PARSE.out.nf_err,
        GENOTYPHI_PARSE.out.nf_log,
        GENOTYPHI_PARSE.out.nf_out,
        GENOTYPHI_PARSE.out.nf_run,
        GENOTYPHI_PARSE.out.nf_sh,
        GENOTYPHI_PARSE.out.nf_trace,
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
