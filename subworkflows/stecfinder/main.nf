//
// stecfinder - Serotype of Shigatoxin producing E. coli using Illumina reads or assemblies
//
nextflow.preview.types = true

include { STECFINDER as STECFINDER_MODULE } from '../../modules/stecfinder/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow STECFINDER {
    take:
    seqs: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    STECFINDER_MODULE(seqs)
    CSVTK_CONCAT(gather(STECFINDER_MODULE.out.tsv, 'stecfinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = STECFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        STECFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
