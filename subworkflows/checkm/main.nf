//
// checkm - Assess the assembly quality of your samples
//
nextflow.preview.types = true

include { CHECKM_LINEAGEWF } from '../../modules/checkm/lineagewf/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CHECKM {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    CHECKM_LINEAGEWF(fasta)

    // Merge results
    ch_merge_checkm = CHECKM_LINEAGEWF.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'checkm'], tsv] }
    CSVTK_CONCAT(ch_merge_checkm, 'tsv', 'tsv')

    emit:
    // Individual outputs
    report = CHECKM_LINEAGEWF.out.tsv
    merged_reports = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = CHECKM_LINEAGEWF.out.tsv.mix(
        CHECKM_LINEAGEWF.out.supplemental,
        CSVTK_CONCAT.out.csv
    )
    logs = CHECKM_LINEAGEWF.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CHECKM_LINEAGEWF.out.nf_begin.mix(
        CHECKM_LINEAGEWF.out.nf_err,
        CHECKM_LINEAGEWF.out.nf_log,
        CHECKM_LINEAGEWF.out.nf_out,
        CHECKM_LINEAGEWF.out.nf_run,
        CHECKM_LINEAGEWF.out.nf_sh,
        CHECKM_LINEAGEWF.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = CHECKM_LINEAGEWF.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
