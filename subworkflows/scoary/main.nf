//
// scoary - GWAS analysis using pangenome outputs
//
nextflow.preview.types = true

include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCOARY {
    take:
    csv: Channel<Tuple<Map, Set<Path>>>
    traits: Path?

    main:    
    SCOARY_MODULE(csv, traits)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.csv
    logs: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.versions
}
