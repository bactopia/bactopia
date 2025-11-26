//
// genotyphi - Salmonella Typhi genotyping with Mykrobe outputs
//
nextflow.preview.types = true

include { MYKROBE_PREDICT } from '../../modules/mykrobe/predict/main'
include { GENOTYPHI_PARSE } from '../../modules/genotyphi/parse/main'
include { CSVTK_CONCAT    } from '../../modules/csvtk/concat/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow GENOTYPHI {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    MYKROBE_PREDICT(reads, "typhi")
    GENOTYPHI_PARSE(MYKROBE_PREDICT.out.json)
    CSVTK_CONCAT(gather(GENOTYPHI_PARSE.out.tsv, 'genotyphi'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = GENOTYPHI_PARSE.out.tsv
    csv: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.csv
    json: Channel<Tuple<Map, Path>> = MYKROBE_PREDICT.out.json
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GENOTYPHI_PARSE.out.tsv,
        MYKROBE_PREDICT.out.csv,
        MYKROBE_PREDICT.out.json,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.logs,
        GENOTYPHI_PARSE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GENOTYPHI_PARSE.out.nf_logs,
        MYKROBE_PREDICT.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MYKROBE_PREDICT.out.versions,
        GENOTYPHI_PARSE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
