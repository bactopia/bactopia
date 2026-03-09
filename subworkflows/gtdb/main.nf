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
 * @input tuple(meta, assembly)
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
 * @output sample_outputs    Per-sample record outputs from CLASSIFY
 * @output run_outputs     Combined GTDB-Tk results across all samples as a record
 */
nextflow.preview.types = true

include { GTDBTK_DOWNLOAD as DOWNLOAD   } from '../../modules/gtdbtk/download/main'
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gather                        } from 'plugin/nf-bactopia'

workflow GTDB {
    take:
    assembly: Channel<Tuple<Map, Path>>
    database: Path
    download_gtdb: Boolean
    save_as_tarball: Boolean

    main:
    if (download_gtdb) {
        // Force CLASSIFY to wait
        DOWNLOAD()

        if (save_as_tarball) {
            CLASSIFY(assembly, DOWNLOAD.out.db_tarball)
        } else {
            CLASSIFY(assembly, DOWNLOAD.out.db)
        }
    } else {
        CLASSIFY(assembly, database)
    }
    CSVTK_CONCAT(gather(CLASSIFY.out, 'gtdb', field: 'bac_tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = CLASSIFY.out
    run_outputs = CSVTK_CONCAT.out
}
