//
// defensefinder - Systematic search of all known anti-phage systems
//
nextflow.preview.types = true

include { DEFENSEFINDER_UPDATE           } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN              } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths                   } from 'plugin/nf-bactopia'
include { gather                         } from 'plugin/nf-bactopia'

workflow DEFENSEFINDER {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]

    main:
    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(fasta, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    GENES_CONCAT(gather(DEFENSEFINDER_RUN.out.genes_tsv, 'defensefinder-genes'), 'tsv', 'tsv')
    HMMER_CONCAT(gather(DEFENSEFINDER_RUN.out.hmmer_tsv, 'defensefinder-hmmer'), 'tsv', 'tsv')
    SYSTEMS_CONCAT(gather(DEFENSEFINDER_RUN.out.systems_tsv, 'defensefinder-systems'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    genes_tsv: Channel<Tuple<Map, Path>> = DEFENSEFINDER_RUN.out.genes_tsv
    merged_genes_tsv: Channel<Tuple<Map, Path>> = GENES_CONCAT.out.csv
    hmmer_tsv: Channel<Tuple<Map, Path>> = DEFENSEFINDER_RUN.out.hmmer_tsv
    merged_hmmer_tsv: Channel<Tuple<Map, Path>> = HMMER_CONCAT.out.csv
    systems_tsv: Channel<Tuple<Map, Path>> = DEFENSEFINDER_RUN.out.systems_tsv
    merged_systems_tsv: Channel<Tuple<Map, Path>> = SYSTEMS_CONCAT.out.csv
    proteins: Channel<Tuple<Map, Path>> = DEFENSEFINDER_RUN.out.proteins
    proteins_index: Channel<Tuple<Map, Path>> = DEFENSEFINDER_RUN.out.proteins_index
    macsydata_raw: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.macsydata_raw

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.genes_tsv,
        GENES_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.hmmer_tsv,
        HMMER_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.systems_tsv,
        SYSTEMS_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.proteins,
        DEFENSEFINDER_RUN.out.proteins_index,
        DEFENSEFINDER_RUN.out.macsydata_raw
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.logs,
        GENES_CONCAT.out.logs,
        HMMER_CONCAT.out.logs,
        SYSTEMS_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.nf_logs,
        GENES_CONCAT.out.nf_logs,
        HMMER_CONCAT.out.nf_logs,
        SYSTEMS_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.versions,
        GENES_CONCAT.out.versions,
        HMMER_CONCAT.out.versions,
        SYSTEMS_CONCAT.out.versions
    ])
}
