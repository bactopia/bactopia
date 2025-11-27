//
// kleborate - Screening Klebsiella genome assemblies for MLST, sub-species, and other Klebsiella related genes of interest
//
nextflow.preview.types = true

include { KLEBORATE as KLEBORATE_MODULE } from '../../modules/kleborate/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow KLEBORATE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    KLEBORATE_MODULE(fasta)
    CSVTK_CONCAT(gather(KLEBORATE_MODULE.out.txt, 'kleborate', 'txt'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = KLEBORATE_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        KLEBORATE_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        KLEBORATE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
