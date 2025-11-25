//
// hicap - Identify cap locus serotype and structure in your Haemophilus influenzae assemblies
//
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow HICAP {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    database_dir
    model_fp

    main:
    HICAP_MODULE(fasta, database_dir, model_fp)

    // Merge results
    ch_merge_hicap = HICAP_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'hicap'], tsv]}
    CSVTK_CONCAT(ch_merge_hicap, 'tsv', 'tsv')

    emit:
    // Individual outputs
    gbk = HICAP_MODULE.out.gbk
    svg = HICAP_MODULE.out.svg
    tsv = HICAP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = HICAP_MODULE.out.gbk.mix(
        HICAP_MODULE.out.svg,
        HICAP_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    )
    logs = HICAP_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
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
    versions = HICAP_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
