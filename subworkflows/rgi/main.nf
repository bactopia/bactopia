//
// rgi - Predict antibiotic resistance from assemblies
//
include { RGI_MAIN } from '../../modules/rgi/main/main'
include { RGI_HEATMAP } from '../../modules/rgi/heatmap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow RGI {
    take:
    fasta // channel: [ val(meta), [ reads ] ]

    main:
    RGI_MAIN(fasta)

    // Merge TSVs
    RGI_MAIN.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'rgi'], tsv]}.set{ ch_merge_rgi }
    CSVTK_CONCAT(ch_merge_rgi, 'tsv', 'tsv')

    // Create Heatmap
    RGI_MAIN.out.json.collect{_meta, json -> json}.map{ json -> [[id:'rgi'], json]}.set{ ch_merge_json }
    RGI_HEATMAP(ch_merge_json)

    emit:
    // Individual outputs
    tsv = RGI_MAIN.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    json = RGI_MAIN.out.json
    heatmap = RGI_HEATMAP.out.heatmap

    // Generic aggregate outputs
    results = RGI_MAIN.out.tsv.mix(
        RGI_HEATMAP.out.heatmap,
        CSVTK_CONCAT.out.csv,
        RGI_MAIN.out.json
    )
    logs = RGI_MAIN.out.logs.mix(
        RGI_HEATMAP.out.logs,
        CSVTK_CONCAT.out.logs
    )
    nf_logs = RGI_MAIN.out.nf_begin.mix(
        RGI_MAIN.out.nf_err,
        RGI_MAIN.out.nf_log,
        RGI_MAIN.out.nf_out,
        RGI_MAIN.out.nf_run,
        RGI_MAIN.out.nf_sh,
        RGI_MAIN.out.nf_trace,
        RGI_HEATMAP.out.nf_begin,
        RGI_HEATMAP.out.nf_err,
        RGI_HEATMAP.out.nf_log,
        RGI_HEATMAP.out.nf_out,
        RGI_HEATMAP.out.nf_run,
        RGI_HEATMAP.out.nf_sh,
        RGI_HEATMAP.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = RGI_MAIN.out.versions.mix(
        RGI_HEATMAP.out.versions,
        CSVTK_CONCAT.out.versions
    )
}
