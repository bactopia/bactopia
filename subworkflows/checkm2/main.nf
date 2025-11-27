//
// checkm2 - Assess the assembly quality of your samples
//
nextflow.preview.types = true

include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT  } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow CHECKM2 {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    database: Path
    download_checkm2: Boolean

    main:
    if (download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(fasta, CHECKM2_DOWNLOAD.out.db)
    } else {
        CHECKM2_PREDICT(fasta, database)
    }

    // Merge results
    CSVTK_CONCAT(gather(CHECKM2_PREDICT.out.tsv, 'checkm2'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = CHECKM2_PREDICT.out.tsv
    merged_reports: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM2_PREDICT.out.tsv,
        CHECKM2_PREDICT.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM2_PREDICT.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM2_PREDICT.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        CHECKM2_PREDICT.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
