//
// seroba - Serotyping of Streptococcus pneumoniae from sequence reads
//
nextflow.preview.types = true

include { SEROBA_RUN   } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow SEROBA {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    SEROBA_RUN(fasta)
    CSVTK_CONCAT(gather(SEROBA_RUN.out.tsv, 'seroba'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SEROBA_RUN.out.tsv
    txt: Channel<Tuple<Map, Path>> = SEROBA_RUN.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.tsv,
        SEROBA_RUN.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SEROBA_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
