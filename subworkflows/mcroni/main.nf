//
// mcroni - Sequence variation in mcr-1 genes (mobilized colistin resistance)
//
nextflow.preview.types = true

include { MCRONI as MCRONI_MODULE } from '../../modules/mcroni/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow MCRONI {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    MCRONI_MODULE(fasta)
    CSVTK_CONCAT(gather(MCRONI_MODULE.out.tsv, 'mcroni'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MCRONI_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    fa: Channel<Tuple<Map, Path>> = MCRONI_MODULE.out.fa

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        MCRONI_MODULE.out.fa
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
