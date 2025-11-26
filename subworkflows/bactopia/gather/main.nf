//
// gather - Tools to gather all samples in one place
//
nextflow.preview.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT            } from '../../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow GATHER {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    GATHER_MODULE(reads)

    // Merge meta values for each sample
    CSVTK_CONCAT(gather(GATHER_MODULE.out.tsv, 'meta'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = GATHER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    fastq_only: Channel<Tuple<Map, Path>> = GATHER_MODULE.out.fastq_only
    raw_fastq: Channel<Tuple<Map, Path>> = GATHER_MODULE.out.raw_fastq
    error: Channel<Tuple<Map, Path>> = GATHER_MODULE.out.error

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.tsv,
        GATHER_MODULE.out.error,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        GATHER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
