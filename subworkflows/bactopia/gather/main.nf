//
// gather - Tools to gather all samples in one place
//
nextflow.preview.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT } from '../../../modules/csvtk/concat/main'

workflow GATHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    GATHER_MODULE(reads)

    // Merge meta values for each sample
    ch_merge_stats = GATHER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'meta'], tsv]}
    CSVTK_CONCAT(ch_merge_stats, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = GATHER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    fastq_only = GATHER_MODULE.out.fastq_only
    raw_fastq = GATHER_MODULE.out.raw_fastq
    error = GATHER_MODULE.out.error

    // Generic aggregate outputs
    results = GATHER_MODULE.out.tsv.mix(
        GATHER_MODULE.out.error,
        CSVTK_CONCAT.out.csv
    )
    logs = GATHER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        GATHER_MODULE.out.nf_begin,
        GATHER_MODULE.out.nf_err,
        GATHER_MODULE.out.nf_log,
        GATHER_MODULE.out.nf_out,
        GATHER_MODULE.out.nf_run,
        GATHER_MODULE.out.nf_sh,
        GATHER_MODULE.out.nf_trace
    )
    versions = GATHER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
