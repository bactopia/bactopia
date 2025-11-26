//
// gamma - A tool for identification, classification, and annotation of translated gene matches
//
nextflow.preview.types = true

include { GAMMA as GAMMA_MODULE } from '../../modules/gamma/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow GAMMA {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]
    db: Channel<Tuple<Map, Path>>

    main:
    GAMMA_MODULE(fasta, db)
    CSVTK_CONCAT(gather(GAMMA_MODULE.out.gamma, 'gamma'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gamma: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.gamma
    merged_gamma: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    psl: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.psl
    fasta: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.fasta
    gff: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.gff

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.gamma,
        CSVTK_CONCAT.out.csv,
        GAMMA_MODULE.out.psl,
        GAMMA_MODULE.out.fasta,
        GAMMA_MODULE.out.gff
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
