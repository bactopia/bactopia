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
    ch_nf_logs = Channel.empty()

    MYKROBE_PREDICT(reads, params.mykrobe_species)
    ch_versions = ch_versions.mix(MYKROBE_PREDICT.out.versions.first())
    ch_logs = ch_logs.mix(MYKROBE_PREDICT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MYKROBE_PREDICT.out.nf_logs)

    // Merge results
    MYKROBE_PREDICT.out.csv.collect{ _meta, csv -> csv }.map{ csv -> [[id:'mykrobe'], csv] }.set{ ch_merge_mykrobe }
    CSVTK_CONCAT(ch_merge_mykrobe, 'csv', 'csv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    csv = MYKROBE_PREDICT.out.csv
    merged_csv = CSVTK_CONCAT.out.csv
    json = MYKROBE_PREDICT.out.json
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
