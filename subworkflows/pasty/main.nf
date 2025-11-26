//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow PASTY {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    PASTY_MODULE(fasta)
    CSVTK_CONCAT(gather(PASTY_MODULE.out.tsv, 'pasty'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = PASTY_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    blast: Channel<Tuple<Map, Path>> = PASTY_MODULE.out.blast

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.tsv,
        PASTY_MODULE.out.blast,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PASTY_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
