//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    PASTY_MODULE(fasta)
    ch_versions = ch_versions.mix(PASTY_MODULE.out.versions)
    ch_logs = ch_logs.mix(PASTY_MODULE.out.logs)

    // Merge results
    PASTY_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'pasty'], tsv]}.set{ ch_merge_pasty }
    CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PASTY_MODULE.out.blast
    logs = ch_logs
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        PASTY_MODULE.out.nf_begin,
        PASTY_MODULE.out.nf_err,
        PASTY_MODULE.out.nf_log,
        PASTY_MODULE.out.nf_out,
        PASTY_MODULE.out.nf_run,
        PASTY_MODULE.out.nf_sh,
        PASTY_MODULE.out.nf_trace
    )
    versions = ch_versions
}
