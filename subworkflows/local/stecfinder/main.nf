//
// stecfinder - Serotype of Shigatoxin producing E. coli using Illumina reads or assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'stecfinder')
options.args = [
    params.stecfinder_use_reads  ? "-r" : "",
    params.stecfinder_hits ? "--hits" : "",
    "--cutoff ${params.stecfinder_cutoff}",
    "--length ${params.stecfinder_length}",
    "--ipaH_length ${params.stecfinder_ipah_length}",
    "--ipaH_depth ${params.stecfinder_ipah_depth}",
    "--stx_length ${params.stecfinder_stx_length}",
    "--stx_depth ${params.stecfinder_stx_depth}",
    "--o_length ${params.stecfinder_o_length}",
    "--o_depth ${params.stecfinder_o_depth}",
    "--h_length ${params.stecfinder_h_length}",
    "--h_depth ${params.stecfinder_h_depth}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { STECFINDER as STECFINDER_MODULE } from '../../../modules/nf-core/stecfinder/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'stecfinder'] )

workflow STECFINDER {
    take:
    seqs // channel: [ val(meta), [ seqs ] ]

    main:
    ch_versions = Channel.empty()

    STECFINDER_MODULE(seqs)
    ch_versions = ch_versions.mix(STECFINDER_MODULE.out.versions.first())

    // Merge results
    STECFINDER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'stecfinder'], tsv]}.set{ ch_merge_stecfinder }
    CSVTK_CONCAT(ch_merge_stecfinder, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = STECFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
