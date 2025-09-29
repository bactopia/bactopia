//
// seqsero2 - Salmonella serotype prediction from reads or assemblies
//
include { SEQSERO2 as SEQSERO2_MODULE } from '../../modules/seqsero2/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SEQSERO2 {
    take:
    seqs // channel: [ val(meta), [ fastqs or assemblies ] ]

    main:
    SEQSERO2_MODULE(seqs)

    // Merge results
    SEQSERO2_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'seqsero2'], tsv]}.set{ ch_merge_seqsero2 }
    CSVTK_CONCAT(ch_merge_seqsero2, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SEQSERO2_MODULE.out.tsv
    txt = SEQSERO2_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = SEQSERO2_MODULE.out.tsv.mix(
        SEQSERO2_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    )
    logs = SEQSERO2_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SEQSERO2_MODULE.out.nf_begin.mix(
        SEQSERO2_MODULE.out.nf_err,
        SEQSERO2_MODULE.out.nf_log,
        SEQSERO2_MODULE.out.nf_out,
        SEQSERO2_MODULE.out.nf_run,
        SEQSERO2_MODULE.out.nf_sh,
        SEQSERO2_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SEQSERO2_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
