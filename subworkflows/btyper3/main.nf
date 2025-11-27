//
// btyper3 - Taxonomic classification of Bacillus cereus group isolates
//
nextflow.preview.types = true

include { BTYPER3 as BTYPER3_MODULE } from '../../modules/btyper3/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow BTYPER3 {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    BTYPER3_MODULE(fasta)
    CSVTK_CONCAT(gather(BTYPER3_MODULE.out.tsv, 'btyper3'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = BTYPER3_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.tsv,
        BTYPER3_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
