//
// busco - Assembly completeness based on evolutionarily informed expectations
// 
nextflow.preview.types = true

include { BUSCO as BUSCO_MODULE } from '../../modules/busco/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow BUSCO {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    busco_lineage: String

    main:
    BUSCO_MODULE(fasta)
    CSVTK_CONCAT(gather(BUSCO_MODULE.out.tsv, "busco-${busco_lineage}"), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = BUSCO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.tsv,
        BUSCO_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
