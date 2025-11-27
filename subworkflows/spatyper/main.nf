//
// SpaTyper - Computational method for finding spa types in Staphylococcus aureus
//
nextflow.preview.types = true

include { SPATYPER as SPATYPER_MODULE } from '../../modules/spatyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SPATYPER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    repeats: Path?
    repeat_order: Path?

    main:
    SPATYPER_MODULE(fasta, repeats, repeat_order)
    CSVTK_CONCAT(gather(SPATYPER_MODULE.out.tsv, 'spatyper'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Path>> = SPATYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate output
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
