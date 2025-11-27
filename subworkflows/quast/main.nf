//
// quast - A module for assessing the quality of assembled contigs
//
nextflow.preview.types = true

include { QUAST as QUAST_MODULE } from '../../modules/quast/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow QUAST {
    take:
    fasta: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    QUAST_MODULE(fasta)
    CSVTK_CONCAT(gather(QUAST_MODULE.out.tsv, 'quast'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = QUAST_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.supplemental,
        QUAST_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        QUAST_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
