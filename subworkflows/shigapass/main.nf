//
// shigapass - Predict Shigella serotypes and differentiate Shigella, EIEC and non-Shigella/EIEC
//
nextflow.preview.types = true

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../modules/shigapass/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow SHIGAPASS {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    SHIGAPASS_MODULE(fasta)

    // Merge results
    ch_merge_shigapass = SHIGAPASS_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'shigapass'], tsv]}
    CSVTK_CONCAT(ch_merge_shigapass, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = SHIGAPASS_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    flex_tsv = SHIGAPASS_MODULE.out.flex_tsv

    // Generic aggregate outputs
    results = SHIGAPASS_MODULE.out.tsv.mix(
        SHIGAPASS_MODULE.out.flex_tsv,
        CSVTK_CONCAT.out.csv
    )
    logs = SHIGAPASS_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = SHIGAPASS_MODULE.out.nf_begin.mix(
        SHIGAPASS_MODULE.out.nf_err,
        SHIGAPASS_MODULE.out.nf_log,
        SHIGAPASS_MODULE.out.nf_out,
        SHIGAPASS_MODULE.out.nf_run,
        SHIGAPASS_MODULE.out.nf_sh,
        SHIGAPASS_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = SHIGAPASS_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
