//
// sccmec - A tool for typing SCCmec cassettes in assemblies
//
include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SCCMEC {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SCCMEC_MODULE(fasta)
    ch_versions = ch_versions.mix(SCCMEC_MODULE.out.versions.first())

    // Merge results
    SCCMEC_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'sccmec'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SCCMEC_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    targets = SCCMEC_MODULE.out.targets
    target_details = SCCMEC_MODULE.out.target_details
    regions = SCCMEC_MODULE.out.regions
    regions_details = SCCMEC_MODULE.out.regions_details
    logs = SCCMEC_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SCCMEC_MODULE.out.nf_begin.mix(
        SCCMEC_MODULE.out.nf_err,
        SCCMEC_MODULE.out.nf_log,
        SCCMEC_MODULE.out.nf_out,
        SCCMEC_MODULE.out.nf_run,
        SCCMEC_MODULE.out.nf_sh,
        SCCMEC_MODULE.out.nf_trace,
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
