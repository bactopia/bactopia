//
// blastn - Search against nucleotide BLAST databases using nucleotide queries
//
nextflow.preview.types = true

include { BLAST_BLASTN as BLASTN_MODULE } from '../../modules/blast/blastn/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTN {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    BLASTN_MODULE(fasta, query)
    CSVTK_CONCAT(gather(BLASTN_MODULE.out.tsv, 'blastn'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = BLASTN_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTN_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
