//
// mlst - Automatic MLST calling from assembled contigs
//
nextflow.preview.types = true

include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT        } from '../../modules/csvtk/concat/main'
include { flattenPaths        } from 'plugin/nf-bactopia'
include { gather              } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    db: Path

    main:
    MLST_MODULE(fasta, db)
    CSVTK_CONCAT(gather(MLST_MODULE.out.tsv, 'mlst'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MLST_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        MLST_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MLST_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
