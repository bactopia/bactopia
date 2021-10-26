//
// hicap - Identify cap locus serotype and structure in your Haemophilus influenzae assemblies
//
hicap_args = [
    params.hicap_debug ? "--debug" : "",
    params.full_sequence  ? "--full_sequence " : "",
    "--gene_coverage ${params.gene_coverage}",
    "--gene_identity ${params.gene_identity}",
    "--broken_gene_length ${params.broken_gene_length}",
    "--broken_gene_identity ${params.broken_gene_identity}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE_DIR = params.database_dir ? file(params.database_dir) : []
MODEL_FP = params.model_fp ? file(params.model_fp) : []

include { HICAP as HICAP_MODULE } from '../../../modules/nf-core/modules/hicap/main' addParams( options: [ args: "${hicap_args}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow HICAP {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    HICAP_MODULE(fasta, DATABASE_DIR, MODEL_FP)
    ch_versions = ch_versions.mix(HICAP_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        HICAP_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'hicap'], tsv]}.set{ ch_merge_hicap }
        CSVTK_CONCAT(ch_merge_hicap, 'tsv', 'tsv')
    }

    emit:
    gbk = HICAP_MODULE.out.gbk
    svg = HICAP_MODULE.out.svg
    tsv = HICAP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
