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
    ch_nf_logs = Channel.empty()

    DATABASE_DIR = params.database_dir ? file(params.database_dir) : []
    MODEL_FP = params.model_fp ? file(params.model_fp) : []

    HICAP_MODULE(fasta, DATABASE_DIR, MODEL_FP)
    ch_versions = ch_versions.mix(HICAP_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(HICAP_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(HICAP_MODULE.out.nf_logs)

    // Merge results
    HICAP_MODULE.out.tsv.collect{ _meta, tsv -> tsv }.map{ tsv -> [[id:'hicap'], tsv] }.set{ ch_merge_hicap }
    CSVTK_CONCAT(ch_merge_hicap, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    gbk = HICAP_MODULE.out.gbk
    svg = HICAP_MODULE.out.svg
    tsv = HICAP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
