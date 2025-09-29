//
// midas - Estimate species abundances from FASTQ files
//
include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MIDAS {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    database

    main:
    MIDAS_SPECIES(reads, database)
    
    // Merge results
    MIDAS_SPECIES.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'midas'], tsv]}.set{ ch_merge_midas }
    CSVTK_CONCAT(ch_merge_midas, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = MIDAS_SPECIES.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    abundances = MIDAS_SPECIES.out.abundances

    // Generic aggregate outputs
    results = MIDAS_SPECIES.out.tsv.mix(
        MIDAS_SPECIES.out.abundances,
        CSVTK_CONCAT.out.csv
    )
    logs = MIDAS_SPECIES.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = MIDAS_SPECIES.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
