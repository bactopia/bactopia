//
// hpsuissero - Serotype prediction of Haemophilus parasuis assemblies
//
nextflow.preview.types = true

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow HPSUISSERO {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    HPSUISSERO_MODULE(fasta)
    CSVTK_CONCAT(gather(HPSUISSERO_MODULE.out.tsv, 'hpsuissero'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = HPSUISSERO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        HPSUISSERO_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        HPSUISSERO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
