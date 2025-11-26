//
// hicap - Identify cap locus serotype and structure in your Haemophilus influenzae assemblies
//
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow HICAP {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]
    database_dir: Channel<Tuple<Map, Path>>
    model_fp: Channel<Tuple<Map, Path>>

    main:
    HICAP_MODULE(fasta, database_dir, model_fp)
    CSVTK_CONCAT(gather(HICAP_MODULE.out.tsv, 'hicap'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gbk: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.gbk
    svg: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.svg
    tsv: Channel<Tuple<Map, Path>> = HICAP_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.gbk,
        HICAP_MODULE.out.svg,
        HICAP_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        HICAP_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
