//
// phispy - Predict prophages in bacterial genomes
//
include { PHISPY as PHISPY_MODULE } from '../../modules/phispy/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PHISPY {
    take:
    gbk // channel: [ val(meta), [ gbk ] ]

    main:
    ch_versions = Channel.empty()

    PHISPY_MODULE(gbk)
    ch_versions = ch_versions.mix(PHISPY_MODULE.out.versions)

    // Merge results
    PHISPY_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'phispy'], tsv]}.set{ ch_merge_phispy }
    CSVTK_CONCAT(ch_merge_phispy, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = PHISPY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = PHISPY_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = PHISPY_MODULE.out.nf_begin.mix(
        PHISPY_MODULE.out.nf_err,
        PHISPY_MODULE.out.nf_log,
        PHISPY_MODULE.out.nf_out,
        PHISPY_MODULE.out.nf_run,
        PHISPY_MODULE.out.nf_sh,
        PHISPY_MODULE.out.nf_trace,
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
