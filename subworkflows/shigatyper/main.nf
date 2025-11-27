//
// shigatyper - Shigella serotype from Illumina or Oxford Nanopore reads
//
nextflow.preview.types = true

include { SHIGATYPER as SHIGATYPER_MODULE } from '../../modules/shigatyper/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow SHIGATYPER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    SHIGATYPER_MODULE(reads)
    CSVTK_CONCAT(gather(SHIGATYPER_MODULE.out.tsv, 'shigatyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SHIGATYPER_MODULE.out.tsv
    hits: Channel<Tuple<Map, Path>> = SHIGATYPER_MODULE.out.hits
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.tsv,
        SHIGATYPER_MODULE.out.hits,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
