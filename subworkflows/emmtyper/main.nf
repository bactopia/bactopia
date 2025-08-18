//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow EMMTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    EMMTYPER_MODULE(fasta, BLASTDB)
    ch_versions = ch_versions.mix(EMMTYPER_MODULE.out.versions.first())
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(EMMTYPER_MODULE.out.logs)

    // Merge results
    EMMTYPER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'emmtyper'], tsv]}.set{ ch_merge_emmtyper }
    CSVTK_CONCAT(ch_merge_emmtyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = EMMTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    BLASTDB = params.emmtyper_blastdb ? file(params.emmtyper_blastdb) : []
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        EMMTYPER_MODULE.out.nf_begin,
        EMMTYPER_MODULE.out.nf_err,
        EMMTYPER_MODULE.out.nf_log,
        EMMTYPER_MODULE.out.nf_out,
        EMMTYPER_MODULE.out.nf_run,
        EMMTYPER_MODULE.out.nf_sh,
        EMMTYPER_MODULE.out.nf_trace
    )
    versions = ch_versions
}
