//
// pirate - Pangenome toolbox for bacterial genomes
//
nextflow.preview.types = true

include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PIRATE {
    take:
    gff: Channel<Tuple<Map, Path>>

    main:
    PIRATE_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.aln
    csv: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PIRATE_MODULE.out.supplemental,
        PIRATE_MODULE.out.aln,
        PIRATE_MODULE.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.versions
}
