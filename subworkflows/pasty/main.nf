//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    PASTY_MODULE(fasta)

    // Merge results
    ch_merge_pasty = PASTY_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'pasty'], tsv]}
    CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = channel.of(CSVTK_CONCAT.out.csv)
    blast = PASTY_MODULE.out.blast

    // Generic aggregate outputs
    results = PASTY_MODULE.out.tsv.mix(
        PASTY_MODULE.out.blast,
        CSVTK_CONCAT.out.csv
    )
    logs = PASTY_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = PASTY_MODULE.out.nf_begin.mix(
        PASTY_MODULE.out.nf_err,
        PASTY_MODULE.out.nf_log,
        PASTY_MODULE.out.nf_out,
        PASTY_MODULE.out.nf_run,
        PASTY_MODULE.out.nf_sh,
        PASTY_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
    )
    versions = PASTY_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
