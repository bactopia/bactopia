//
// seroba - Serotyping of Streptococcus pneumoniae from sequence reads
//
include { SEROBA_RUN } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SEROBA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SEROBA_RUN(fasta)
    ch_versions = ch_versions.mix(SEROBA_RUN.out.versions.first())

    // Merge results
    SEROBA_RUN.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'seroba'], tsv]}.set{ ch_merge_seroba }
    CSVTK_CONCAT(ch_merge_seroba, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SEROBA_RUN.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = SEROBA_RUN.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SEROBA_RUN.out.nf_begin.mix(
        SEROBA_RUN.out.nf_err,
        SEROBA_RUN.out.nf_log,
        SEROBA_RUN.out.nf_out,
        SEROBA_RUN.out.nf_run,
        SEROBA_RUN.out.nf_sh,
        SEROBA_RUN.out.nf_trace,
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
