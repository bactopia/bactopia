//
// btyper3 - Taxonomic classification of Bacillus cereus group isolates
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'btyper3')
options.args = [
    params.bt_opts ? "${params.bt_opts}" : "",
    "--virulence_identity ${params.bt_virulence_identity}",
    "--virulence_coverage ${params.bt_virulence_coverage}",
    "--bt_identity ${params.bt_identity}",
    "--bt_coverage ${params.bt_coverage}",
    "--bt_overlap ${params.bt_overlap}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { BTYPER3 as BTYPER3_MODULE } from '../../../modules/nf-core/btyper3/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'btyper3-concat', process_name: params.merge_folder] )

workflow BTYPER3 {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    BTYPER3_MODULE(fasta)
    ch_versions = ch_versions.mix(BTYPER3_MODULE.out.versions.first())

    // Merge results
    BTYPER3_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'btyper3'], tsv]}.set{ ch_merge_btyper3 }
    CSVTK_CONCAT(ch_merge_btyper3, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = BTYPER3_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
