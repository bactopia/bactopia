//
// ssuissero - Serotype prediction of Streptococcus suis assemblies
//
nextflow.preview.types = true

include { SSUISSERO as SSUISSERO_MODULE } from '../../modules/ssuissero/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow SSUISSERO {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    SSUISSERO_MODULE(fasta)
    CSVTK_CONCAT(gather(SSUISSERO_MODULE.out.tsv, 'ssuissero'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SSUISSERO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
