//
// seqsero2 - Salmonella serotype prediction from reads or assemblies
//
nextflow.preview.types = true

include { SEQSERO2 as SEQSERO2_MODULE } from '../../modules/seqsero2/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SEQSERO2 {
    take:
    seqs: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fastqs or assemblies ] ]

    main:
    SEQSERO2_MODULE(seqs)
    CSVTK_CONCAT(gather(SEQSERO2_MODULE.out.tsv, 'seqsero2'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SEQSERO2_MODULE.out.tsv
    txt: Channel<Tuple<Map, Path>> = SEQSERO2_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.tsv,
        SEQSERO2_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
