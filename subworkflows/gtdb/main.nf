//
// gtdb - Identify marker genes and assign taxonomic classifications
//
nextflow.preview.types = true

include { GTDBTK_DOWNLOAD as DOWNLOAD   } from '../../modules/gtdbtk/download/main'
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../modules/gtdbtk/classifywf/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow GTDB {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    database: Path
    download_gtdb: Boolean
    save_as_tarball: Boolean

    main:
    if (download_gtdb) {
        // Force CLASSIFY to wait
        DOWNLOAD()

        if (save_as_tarball) {
            CLASSIFY(fasta, DOWNLOAD.out.db_tarball)
        } else {
            CLASSIFY(fasta, DOWNLOAD.out.db)
        }
    } else {
        CLASSIFY(fasta, database)
    }
    CSVTK_CONCAT(gather(CLASSIFY.out.tsv, 'gtdb'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = CLASSIFY.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

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
