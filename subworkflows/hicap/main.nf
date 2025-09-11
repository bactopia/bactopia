//
// hicap - Identify cap locus serotype and structure in your Haemophilus influenzae assemblies
//
include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow HICAP {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    DATABASE_DIR = params.database_dir ? file(params.database_dir) : []
    MODEL_FP = params.model_fp ? file(params.model_fp) : []

    HICAP_MODULE(fasta, DATABASE_DIR, MODEL_FP)
    ch_versions = ch_versions.mix(HICAP_MODULE.out.versions)
    ch_logs = ch_logs.mix(HICAP_MODULE.out.logs)

    // Aggregate results
    CSVTK_CONCAT(
        HICAP_MODULE.out.tsv.collect{it[1]}.map{ tsv -> [[id: 'hicap'], tsv]},
        'tsv',
        'tsv'
    )
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    gbk = HICAP_MODULE.out.gbk
    svg = HICAP_MODULE.out.svg
    tsv = HICAP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = HICAP_MODULE.out.nf_begin.mix(
        HICAP_MODULE.out.nf_err,
        HICAP_MODULE.out.nf_log,
        HICAP_MODULE.out.nf_out,
        HICAP_MODULE.out.nf_run,
        HICAP_MODULE.out.nf_sh,
        HICAP_MODULE.out.nf_trace,
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
