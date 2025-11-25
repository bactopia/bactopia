//
// meningotype - Serotyping of Neisseria meningitidis
//
nextflow.preview.types = true

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    MENINGOTYPE_MODULE(fasta)

    // Merge results
    ch_merge_meningotype = MENINGOTYPE_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'meningotype'], tsv]}
    CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = MENINGOTYPE_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = MENINGOTYPE_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        MENINGOTYPE_MODULE.out.nf_begin,
        MENINGOTYPE_MODULE.out.nf_err,
        MENINGOTYPE_MODULE.out.nf_log,
        MENINGOTYPE_MODULE.out.nf_out,
        MENINGOTYPE_MODULE.out.nf_run,
        MENINGOTYPE_MODULE.out.nf_sh,
        MENINGOTYPE_MODULE.out.nf_trace
    )
    versions = MENINGOTYPE_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
