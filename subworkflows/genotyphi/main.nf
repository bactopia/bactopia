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
    ch_nf_logs = Channel.empty()

    MYKROBE_PREDICT(reads, "typhi")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())
    ch_logs = ch_logs.mix(MYKROBE_PREDICT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MYKROBE_PREDICT.out.nf_logs)

    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    ch_versions = ch_versions.mix(GENOTYPHI_PARSE.out.versions.first())
    ch_logs = ch_logs.mix(GENOTYPHI_PARSE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GENOTYPHI_PARSE.out.nf_logs)

    // Merge results
    GENOTYPHI_PARSE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'genotyphi'], tsv] }.set{ ch_merge_genotyphi }
    CSVTK_CONCAT(ch_merge_genotyphi, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = GENOTYPHI_PARSE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    csv = MYKROBE_PREDICT.out.csv
    json = MYKROBE_PREDICT.out.json
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
