//
// meningotype - Serotyping of Neisseria meningitidis
//
nextflow.preview.types = true

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../modules/meningotype/main'
include { CSVTK_CONCAT                      } from '../../modules/csvtk/concat/main'
include { flattenPaths                      } from 'plugin/nf-bactopia'
include { gather                            } from 'plugin/nf-bactopia'

workflow MENINGOTYPE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    MENINGOTYPE_MODULE(fasta)
    CSVTK_CONCAT(gather(MENINGOTYPE_MODULE.out.tsv, 'meningotype'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MENINGOTYPE_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        MENINGOTYPE_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MENINGOTYPE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
