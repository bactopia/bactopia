//
// fastani - fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
//
nextflow.preview.types = true

include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow FASTANI {
    take:
    query: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]
    reference: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    FASTANI_MODULE(gather(query, 'query', 'fasta'), gather(reference, 'reference', 'fasta'))
    CSVTK_CONCAT(gather(FASTANI_MODULE.out.tsv, 'fastani'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = FASTANI_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
