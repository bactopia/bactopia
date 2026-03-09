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
 * @output sample_outputs  Per-sample record outputs from CHECKM2_PREDICT
 * @output run_outputs   Combined CheckM2 results across all samples as a record
 */
nextflow.preview.types = true

include { CHECKM2_DOWNLOAD } from '../../modules/checkm2/download/main'
include { CHECKM2_PREDICT  } from '../../modules/checkm2/predict/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gather           } from 'plugin/nf-bactopia'

workflow CHECKM2 {
    take:
    assembly: Channel<Tuple<Map, Path>>
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
    CSVTK_CONCAT(gather(CHECKM2_PREDICT.out, 'checkm2', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = CHECKM2_PREDICT.out
    run_outputs = CSVTK_CONCAT.out
}
