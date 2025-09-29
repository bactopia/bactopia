//
// quast - A module for assessing the quality of assembled contigs
//
include { QUAST as QUAST_MODULE } from '../../modules/quast/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow QUAST {
    take:
    fasta // channel: [ val(meta), [ fasta ], [ meta_files ] ]

    main:
    QUAST_MODULE(fasta)

    // Merge results
    QUAST_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'quast'], tsv]}.set{ ch_merge_quast}
    CSVTK_CONCAT(ch_merge_quast, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = QUAST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = QUAST_MODULE.out.supplemental.mix(
        QUAST_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    )
    logs = QUAST_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = QUAST_MODULE.out.nf_begin.mix(
        QUAST_MODULE.out.nf_err,
        QUAST_MODULE.out.nf_log,
        QUAST_MODULE.out.nf_out,
        QUAST_MODULE.out.nf_run,
        QUAST_MODULE.out.nf_sh,
        QUAST_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = QUAST_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
