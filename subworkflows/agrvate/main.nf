//
// agrvate - Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
//
nextflow.preview.types = true

include { AGRVATE as AGRVATE_MODULE } from '../../modules/agrvate/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow AGRVATE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    AGRVATE_MODULE(fasta)
    CSVTK_CONCAT(gather(AGRVATE_MODULE.out.summary, 'agrvate'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Path>> = AGRVATE_MODULE.out.summary
    supplemental: Channel<Tuple<Map, Path>> = AGRVATE_MODULE.out.supplemental
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate output
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.summary,
        AGRVATE_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
