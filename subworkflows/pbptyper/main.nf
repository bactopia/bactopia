//
// pbptyper - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae
//
include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PBPTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    PBPTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(PBPTYPER_MODULE.out.versions.first()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions))
    ch_logs = ch_logs.mix(PBPTYPER_MODULE.out.logs)

    // Merge results
    PBPTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'pbptyper'], tsv]}.set{ ch_merge_pbptyper }
    CSVTK_CONCAT(ch_merge_pbptyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = PBPTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PBPTYPER_MODULE.out.blast
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        PBPTYPER_MODULE.out.nf_begin,
        PBPTYPER_MODULE.out.nf_err,
        PBPTYPER_MODULE.out.nf_log,
        PBPTYPER_MODULE.out.nf_out,
        PBPTYPER_MODULE.out.nf_run,
        PBPTYPER_MODULE.out.nf_sh,
        PBPTYPER_MODULE.out.nf_trace
    )
    versions = ch_versions
}
