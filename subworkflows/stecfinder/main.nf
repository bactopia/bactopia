//
// stecfinder - Serotype of Shigatoxin producing E. coli using Illumina reads or assemblies
//
include { STECFINDER as STECFINDER_MODULE } from '../../modules/stecfinder/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow STECFINDER {
    take:
    seqs // channel: [ val(meta), [ seqs ] ]

    main:
    ch_versions = Channel.empty()

    STECFINDER_MODULE(seqs)
    ch_versions = ch_versions.mix(STECFINDER_MODULE.out.versions)

    // Merge results
    STECFINDER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'stecfinder'], tsv]}.set{ ch_merge_stecfinder }
    CSVTK_CONCAT(ch_merge_stecfinder, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = STECFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = STECFINDER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = STECFINDER_MODULE.out.nf_begin.mix(
        STECFINDER_MODULE.out.nf_err,
        STECFINDER_MODULE.out.nf_log,
        STECFINDER_MODULE.out.nf_out,
        STECFINDER_MODULE.out.nf_run,
        STECFINDER_MODULE.out.nf_sh,
        STECFINDER_MODULE.out.nf_trace,
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
