//
// tblastn - Search against translated nucleotide BLAST databases using protein queries
//
nextflow.preview.types = true

include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../modules/blast/tblastn/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow TBLASTN {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    TBLASTN_MODULE(fasta, query)
    CSVTK_CONCAT(gather(TBLASTN_MODULE.out.tsv, 'tblastn'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = TBLASTN_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        TBLASTN_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
