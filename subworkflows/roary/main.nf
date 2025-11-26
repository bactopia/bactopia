//
// roary - Rapid large-scale prokaryote pangenome analysis
//
nextflow.preview.types = true

include { ROARY as ROARY_MODULE } from '../../modules/roary/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow ROARY {
    take:
    gff : Channel<Tuple<Map, List<Path>>>

    main:
    ROARY_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Path>> = ROARY_MODULE.out.aln
    csv: Channel<Tuple<Map, Path>> = ROARY_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ROARY_MODULE.out.supplemental,
        ROARY_MODULE.out.aln,
        ROARY_MODULE.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = ROARY_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = ROARY_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = ROARY_MODULE.out.versions
}
