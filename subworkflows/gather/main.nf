//
// gather - Tools to gather all samples in one place
//
include { GATHER as GATHER_MODULE } from '../../modules/bactopia/gather/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow GATHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    // Gather genomes (local, assembly, SRA/ENA)
    GATHER_MODULE(reads)
    ch_versions = ch_versions.mix(GATHER_MODULE.out.versions)
    ch_logs = ch_logs.mix(GATHER_MODULE.out.logs)
    // Merge meta values for each sample
    GATHER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'meta'], tsv]}.set{ ch_merge_stats }
    CSVTK_CONCAT(ch_merge_stats, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:

    emit:
    tsv = GATHER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    fastq_only = GATHER_MODULE.out.fastq_only
    raw_fastq = GATHER_MODULE.out.raw_fastq
    logs = ch_logs
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
    versions = ch_versions
}
