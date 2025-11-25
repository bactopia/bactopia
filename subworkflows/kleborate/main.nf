//
// kleborate - Screening Klebsiella genome assemblies for MLST, sub-species, and other Klebsiella related genes of interest
//
nextflow.preview.types = true

include { KLEBORATE as KLEBORATE_MODULE } from '../../modules/kleborate/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow KLEBORATE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    KLEBORATE_MODULE(fasta)

    // Merge results
    ch_merge_kleborate = KLEBORATE_MODULE.out.txt.collect{_meta, txt -> txt}.map{ txt -> [[id:'kleborate'], txt]}
    CSVTK_CONCAT(ch_merge_kleborate, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = KLEBORATE_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = KLEBORATE_MODULE.out.txt.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = KLEBORATE_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        KLEBORATE_MODULE.out.nf_begin,
        KLEBORATE_MODULE.out.nf_err,
        KLEBORATE_MODULE.out.nf_log,
        KLEBORATE_MODULE.out.nf_out,
        KLEBORATE_MODULE.out.nf_run,
        KLEBORATE_MODULE.out.nf_sh,
        KLEBORATE_MODULE.out.nf_trace
    )
    versions = KLEBORATE_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
