//
// mykrobe - Antimicrobial resistance detection for specific species
//
nextflow.preview.types = true

include { MYKROBE_PREDICT }  from '../../modules/mykrobe/predict/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow MYKROBE {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]
    mykrobe_species: String

    main:
    MYKROBE_PREDICT(reads, mykrobe_species)
    CSVTK_CONCAT(gather(MYKROBE_PREDICT.out.csv, 'mykrobe'), 'csv', 'csv')

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.csv
    json: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.json
    merged_csv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.csv,
        MYKROBE_PREDICT.out.json,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
