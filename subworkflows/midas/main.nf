//
// midas - Estimate species abundances from FASTQ files
//
nextflow.preview.types = true

include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT  } from '../../modules/csvtk/concat/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow MIDAS {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]
    database: Channel<Tuple<Map, Path>>

    main:
    MIDAS_SPECIES(reads, database)
    CSVTK_CONCAT(gather(MIDAS_SPECIES.out.tsv, 'midas'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MIDAS_SPECIES.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    abundances: Channel<Tuple<Map, Path>> = MIDAS_SPECIES.out.abundances

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.tsv,
        MIDAS_SPECIES.out.abundances,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
