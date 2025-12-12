/**
 * Assess metagenome bin completeness using CheckM2.
 *
 * This subworkflow evaluates the quality and completeness of metagenome-assembled genomes
 * (MAGs) using [CheckM2](https://github.com/chklovski/CheckM2). It provides an improved
 * assessment using machine learning models trained on high-quality reference genomes,
 * offering more accurate completeness and contamination estimates. The workflow can either
 * download the required database or use a user-provided database path.
 *
 * @status stable
 * @keywords metagenome, bin, completeness, contamination, mag, quality, machine-learning
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, conditional-logic, database-dependent
 * @citation checkm2
 *
 * @modules csvtk_concat, checkm2_predict, checkm2_download
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Metagenome-assembled genome bins to evaluate. Each tuple contains metadata
 *   about the sample and a set of genome bins in FASTA format.
 *
 * @input database
 * Path to CheckM2 database directory. If download_checkm2 is true, this can be
 * a placeholder as the database will be downloaded automatically.
 *
 * @input download_checkm2
 * Boolean flag to automatically download the CheckM2 database if not available.
 * When true, downloads the required reference database before prediction.
 *
 * @output report        Per-bin CheckM2 quality assessment results in TSV format
 * @output merged_reports Combined CheckM2 results summary across all bins
 * @output results       Aggregated results channel containing all output files
 * @output logs          Aggregated logs channel containing all execution logs
 * @output nf_logs       Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions      Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT  } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow CHECKM2 {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    database: Path
    download_checkm2: Boolean

    main:
    if (download_checkm2) {
        CHECKM2_DOWNLOAD()
        CHECKM2_PREDICT(assembly, CHECKM2_DOWNLOAD.out.db)
    } else {
        CHECKM2_PREDICT(assembly, database)
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
