//
// ectyper - In-silico prediction of Escherichia coli serotype
//
nextflow.preview.types = true

include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow ECTYPER {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    ECTYPER_MODULE(fasta)
    CSVTK_CONCAT(gather(ECTYPER_MODULE.out.tsv, 'ectyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = ECTYPER_MODULE.out.tsv
    txt: Channel<Tuple<Map, Path>> = ECTYPER_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.tsv,
        ECTYPER_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        ECTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ECTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
