//
// assembler - Assembly of Illumina and ONT reads
//
include { ASSEMBLER as ASSEMBLER_MODULE } from '../../modules/bactopia/assembler/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ASSEMBLER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ASSEMBLER_MODULE(reads)

    // Merge results
    ASSEMBLER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'assembly-scan'], tsv]}.set{ ch_merge_stats }
    CSVTK_CONCAT(ch_merge_stats, 'tsv', 'tsv')

    emit:
    // Individual outputs
    fna = ASSEMBLER_MODULE.out.fna
    tsv = ASSEMBLER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = ASSEMBLER_MODULE.out.fna.mix(
        ASSEMBLER_MODULE.out.tsv,
        ASSEMBLER_MODULE.out.error,
        CSVTK_CONCAT.out.csv
    )
    logs = ASSEMBLER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = ASSEMBLER_MODULE.out.nf_begin.mix(
        ASSEMBLER_MODULE.out.nf_err,
        ASSEMBLER_MODULE.out.nf_log,
        ASSEMBLER_MODULE.out.nf_out,
        ASSEMBLER_MODULE.out.nf_run,
        ASSEMBLER_MODULE.out.nf_sh,
        ASSEMBLER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ASSEMBLER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
