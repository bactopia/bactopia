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
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation,conditional-logic,database-dependent
 * @citation checkm2
 *
 * @modules csvtk_concat, checkm2_predict, checkm2_download
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Metagenome-assembled genome bins to evaluate. Each record contains metadata
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
 * @output sample_outputs
 * - `tsv`: A tab-delimited report of quality metrics (Completeness, Contamination)
 * - `supplemental`: Directory containing intermediate protein files and Diamond alignments
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT  } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gatherCsvtk      } from 'plugin/nf-bactopia'

workflow CHECKM2 {
    take:
    assembly: Channel<Record>
    database: Path
    download_checkm2: Boolean

    main:
    ch_checkm2_predict = channel.empty()
    if (download_checkm2) {
        ch_checkm2_download = CHECKM2_DOWNLOAD()
        ch_checkm2_predict = CHECKM2_PREDICT(assembly, ch_checkm2_download.map { r -> r.db })
    } else {
        ch_checkm2_predict = CHECKM2_PREDICT(assembly, database)
    }

    // Merge results
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_checkm2_predict, 'tsv', [name: 'checkm2']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_checkm2_predict
    run_outputs = ch_csvtk_concat
}
