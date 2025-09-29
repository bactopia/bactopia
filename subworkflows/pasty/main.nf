//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    PASTY_MODULE(fasta)

    // Merge results
    PASTY_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'pasty'], tsv]}.set{ ch_merge_pasty }
    CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PASTY_MODULE.out.blast

    // Generic aggregate outputs
    results = PASTY_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv,
        PASTY_MODULE.out.blast
    )
    logs = PASTY_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = PASTY_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
