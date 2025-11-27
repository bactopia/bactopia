//
// shigeifinder - Shigella and EIEC serotyping from assemblies
//
nextflow.preview.types = true

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../modules/shigeifinder/main'
include { CSVTK_CONCAT                        } from '../../modules/csvtk/concat/main'
include { flattenPaths                        } from 'plugin/nf-bactopia'
include { gather                              } from 'plugin/nf-bactopia'

workflow SHIGEIFINDER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    SHIGEIFINDER_MODULE(fasta)
    CSVTK_CONCAT(gather(SHIGEIFINDER_MODULE.out.tsv, 'shigeifinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SHIGEIFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGEIFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
