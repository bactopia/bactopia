//
// shigapass - Predict Shigella serotypes and differentiate Shigella, EIEC and non-Shigella/EIEC
//
nextflow.preview.types = true

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../modules/shigapass/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow SHIGAPASS {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    SHIGAPASS_MODULE(fasta)
    CSVTK_CONCAT(gather(SHIGAPASS_MODULE.out.tsv, 'shigapass'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SHIGAPASS_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    flex_tsv: Channel<Tuple<Map, Path>> = SHIGAPASS_MODULE.out.flex_tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.tsv,
        SHIGAPASS_MODULE.out.flex_tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
