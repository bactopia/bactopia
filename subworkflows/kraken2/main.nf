//
// kraken2 - Taxonomic classification of sequence reads 
//
nextflow.preview.types = true

include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow KRAKEN2 {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    // Individual outputs
    classified: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.classified
    kraken2_report: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.kraken2_report
    unclassified: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.unclassified

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2_MODULE.out.classified,
        KRAKEN2_MODULE.out.kraken2_report,
        KRAKEN2_MODULE.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.versions
}
