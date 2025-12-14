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
 * @output tsv         GTDB-Tk classification results with taxonomic assignments
 * @output merged_tsv  Combined TSV file containing classification results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GTDBTK_DOWNLOAD as DOWNLOAD   } from '../../modules/gtdbtk/download/main'
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
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
    CSVTK_CONCAT(gather(CLASSIFY.out.tsv, 'gtdb'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = CLASSIFY.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        CLASSIFY.out.tsv,
        CLASSIFY.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CLASSIFY.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CLASSIFY.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        CLASSIFY.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
