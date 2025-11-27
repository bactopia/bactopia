//
// sccmec - A tool for typing SCCmec cassettes in assemblies
//
nextflow.preview.types = true

include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCCMEC {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    SCCMEC_MODULE(fasta)
    CSVTK_CONCAT(gather(SCCMEC_MODULE.out.tsv, 'sccmec'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    targets: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.targets
    target_details: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.target_details
    regions: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.regions
    regions_details: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.regions_details

    // Generic aggregate output
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        SCCMEC_MODULE.out.targets,
        SCCMEC_MODULE.out.target_details,
        SCCMEC_MODULE.out.regions,
        SCCMEC_MODULE.out.regions_details
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
