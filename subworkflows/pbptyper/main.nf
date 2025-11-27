//
// pbptyper - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae
//
nextflow.preview.types = true

include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow PBPTYPER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    PBPTYPER_MODULE(fasta)
    CSVTK_CONCAT(gather(PBPTYPER_MODULE.out.tsv, 'pbptyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = PBPTYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    blast: Channel<Tuple<Map, Path>> = PBPTYPER_MODULE.out.blast

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PBPTYPER_MODULE.out.blast
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        PBPTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PBPTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
