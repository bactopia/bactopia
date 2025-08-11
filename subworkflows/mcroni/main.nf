//
// mcroni - Sequence variation in mcr-1 genes (mobilized colistin resistance)
//
include { MCRONI as MCRONI_MODULE } from '../../modules/mcroni/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MCRONI {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MCRONI_MODULE(fasta)
    ch_versions = ch_versions.mix(MCRONI_MODULE.out.versions.first())

    // Merge results
    MCRONI_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'mcroni'], tsv]}.set{ ch_merge_mcroni }
    CSVTK_CONCAT(ch_merge_mcroni, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = MCRONI_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    fa = MCRONI_MODULE.out.fa
    logs = MCRONI_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = MCRONI_MODULE.out.nf_begin.mix(
        MCRONI_MODULE.out.nf_err,
        MCRONI_MODULE.out.nf_log,
        MCRONI_MODULE.out.nf_out,
        MCRONI_MODULE.out.nf_run,
        MCRONI_MODULE.out.nf_sh,
        MCRONI_MODULE.out.nf_trace,
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
