//
// checkm - Assess the assembly quality of your samples
//
nextflow.preview.types = true

include { CHECKM_LINEAGEWF } from '../../modules/checkm/lineagewf/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow CHECKM {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    CHECKM_LINEAGEWF(fasta)
    CSVTK_CONCAT(gather(CHECKM_LINEAGEWF.out.tsv, 'checkm'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = CHECKM_LINEAGEWF.out.tsv
    merged_reports: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.tsv,
        CHECKM_LINEAGEWF.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM_LINEAGEWF.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
