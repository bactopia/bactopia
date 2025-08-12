//
// midas - Estimate species abundances from FASTQ files
//
include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MIDAS {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    DATABASE = params.midas_db ? file(params.midas_db) : []

    MIDAS_SPECIES(reads, DATABASE)
    ch_versions = ch_versions.mix(MIDAS_SPECIES.out.versions)
    ch_logs = ch_logs.mix(MIDAS_SPECIES.out.logs)
    
    // Merge results
    MIDAS_SPECIES.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'midas'], tsv]}.set{ ch_merge_midas }
    CSVTK_CONCAT(ch_merge_midas, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = MIDAS_SPECIES.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    abundances = MIDAS_SPECIES.out.abundances
    logs = ch_logs
    nf_logs = MIDAS_SPECIES.out.nf_begin.mix(
        MIDAS_SPECIES.out.nf_err,
        MIDAS_SPECIES.out.nf_log,
        MIDAS_SPECIES.out.nf_out,
        MIDAS_SPECIES.out.nf_run,
        MIDAS_SPECIES.out.nf_sh,
        MIDAS_SPECIES.out.nf_trace,
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
