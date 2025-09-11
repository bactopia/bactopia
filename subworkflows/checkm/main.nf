//
// checkm - Assess the assembly quality of your samples
//
include { CHECKM_LINEAGEWF } from '../../modules/checkm/lineagewf/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CHECKM {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_merged_checkm = Channel.empty()

    CHECKM_LINEAGEWF(fasta)
    ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
    ch_logs = ch_logs.mix(CHECKM_LINEAGEWF.out.logs)

    CHECKM_LINEAGEWF.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'checkm'], tsv] }.set{ ch_merge_checkm }
    CSVTK_CONCAT(ch_merge_checkm, 'tsv', 'tsv')
    ch_merged_checkm = ch_merged_checkm.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    results = CHECKM_LINEAGEWF.out.results
    report = CHECKM_LINEAGEWF.out.tsv
    merged_reports = ch_merged_checkm
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = CHECKM_LINEAGEWF.out.nf_begin.mix(
        CHECKM_LINEAGEWF.out.nf_err,
        CHECKM_LINEAGEWF.out.nf_log,
        CHECKM_LINEAGEWF.out.nf_out,
        CHECKM_LINEAGEWF.out.nf_run,
        CHECKM_LINEAGEWF.out.nf_sh,
        CHECKM_LINEAGEWF.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
