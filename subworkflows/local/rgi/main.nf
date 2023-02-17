//
// rgi - Predict antibiotic resistance from assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'rgi')
options.args = [
    params.use_diamond ? "--alignment_tool DIAMOND" : "--alignment_tool BLAST",
    params.include_loose ? "--include_loose" : "",
    params.exclude_nudge ? "--exclude_nudge" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

heatmap_args = [
    params.rgi_frequency ? "--frequency" : "",
    params.rgi_category ? "--category ${params.rgi_category}" : "",
    params.rgi_cluster ? "--cluster ${params.rgi_cluster}" : "",
    "--display ${params.rgi_display}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { RGI_MAIN } from '../../../modules/nf-core/rgi/main/main' addParams( options: options )
include { RGI_HEATMAP } from '../../../modules/nf-core/rgi/heatmap/main' addParams( options: [args: heatmap_args, logs_subdir: 'rgi-heatmap', process_name: params.merge_folder] )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'rgi-concat', process_name: params.merge_folder] )

workflow RGI {
    take:
    fasta // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    RGI_MAIN(fasta)
    ch_versions = ch_versions.mix(RGI_MAIN.out.versions.first())

    // Merge TSVs
    RGI_MAIN.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'rgi'], tsv]}.set{ ch_merge_rgi }
    CSVTK_CONCAT(ch_merge_rgi, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    // Create Heatmap
    RGI_MAIN.out.json.collect{meta, json -> json}.map{ json -> [[id:'rgi'], json]}.set{ ch_merge_json }
    RGI_HEATMAP(ch_merge_json)
    ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

    emit:
    tsv = RGI_MAIN.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    json = RGI_MAIN.out.json
    heatmap = RGI_HEATMAP.out.heatmap
    versions = ch_versions // channel: [ versions.yml ]
}
