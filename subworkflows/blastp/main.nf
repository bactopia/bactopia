//
// blastp - Search against protein BLAST databases using protein queries
//
nextflow.preview.types = true

include { BLAST_BLASTP as BLASTP_MODULE } from '../../modules/blast/blastp/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow BLASTP {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    query: Path

    main:
    BLASTP_MODULE(fasta, query)
    CSVTK_CONCAT(gather(BLASTP_MODULE.out.tsv, 'blastp'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = BLASTP_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTP_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTP_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTP_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BLASTP_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
