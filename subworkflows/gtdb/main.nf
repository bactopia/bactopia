/**
 * Taxonomic classification with the Genome Taxonomy Database.
 *
 * This subworkflow assigns objective taxonomic classifications to bacterial and
 * archaeal genomes using [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk), which
 * is based on the Genome Taxonomy Database (GTDB). The workflow can optionally
 * download the GTDB database and supports both unpacked and tarball database formats.
 * It provides taxonomic placement and phylogenetic marker gene identification.
 *
 * @status stable
 * @keywords taxonomy, classification, GTDB, phylogeny, marker genes
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic, resource-download
 * @citation gtdb_tk
 *
 * @modules csvtk_concat, gtdbtk_classifywf, gtdbtk_download
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for taxonomic classification
 *
 * @input database
 * Path to GTDB reference database (optional, will download if not provided)
 *
 * @input download_gtdb
 * Boolean flag to trigger GTDB database download if not provided
 *
 * @input save_as_tarball
 * Boolean flag to use tarball format database when downloading
 *
 * @output sample_outputs
 * - `bac_tsv`: The bacterial classification summary file containing the taxonomic assignment
 * - `ar_tsv`: The archaeal classification summary file containing the taxonomic assignment
 * - `supplemental`: Directory containing the reference tree, alignments, and detailed logs
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { GTDBTK_DOWNLOAD   } from '../../modules/gtdbtk/download/main'
include { GTDBTK_CLASSIFYWF } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT      } from '../../modules/csvtk/concat/main'
include { gatherCsvtk            } from 'plugin/nf-bactopia'

workflow GTDB {
    take:
    assembly: Channel<Record>
    database: Path
    download_gtdb: Boolean
    save_as_tarball: Boolean

    main:
    if (download_gtdb) {
        // Force CLASSIFY to wait
        GTDBTK_DOWNLOAD()

        if (save_as_tarball) {
            GTDBTK_CLASSIFYWF(assembly, GTDBTK_DOWNLOAD.out.map { r -> r.db_tarball })
        } else {
            GTDBTK_CLASSIFYWF(assembly, GTDBTK_DOWNLOAD.out.map { r -> r.db })
        }
    } else {
        GTDBTK_CLASSIFYWF(assembly, database)
    }
    CSVTK_CONCAT(gatherCsvtk(GTDBTK_CLASSIFYWF.out, 'bac_tsv', [name: 'gtdb']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = GTDBTK_CLASSIFYWF.out
    run_outputs = CSVTK_CONCAT.out
}
