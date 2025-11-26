//
// rgi - Predict antibiotic resistance from assemblies
//
nextflow.preview.types = true

include { RGI_MAIN     } from '../../modules/rgi/main/main'
include { RGI_HEATMAP  } from '../../modules/rgi/heatmap/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow RGI {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    RGI_MAIN(fasta)
    CSVTK_CONCAT(gather(RGI_MAIN.out.tsv, 'rgi'), 'tsv', 'tsv')
    RGI_HEATMAP(gather(RGI_MAIN.out.json, 'rgi', 'json'))

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = RGI_MAIN.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    json: Channel<Tuple<Map, Path>> = RGI_MAIN.out.json
    heatmap: Channel<Tuple<Map, Path>> = RGI_HEATMAP.out.heatmap

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.tsv,
        RGI_HEATMAP.out.heatmap,
        CSVTK_CONCAT.out.csv,
        RGI_MAIN.out.json
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.logs,
        RGI_HEATMAP.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.nf_logs,
        RGI_HEATMAP.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        RGI_MAIN.out.versions,
        RGI_HEATMAP.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
