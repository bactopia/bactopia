//
// emmtyper - emm-typing of Streptococcus pyogenes assemblies
//
nextflow.preview.types = true

include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow EMMTYPER {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ fasta ] ]
    blastdb: Channel<Tuple<Map, Path>>

    main:
    EMMTYPER_MODULE(fasta, blastdb)
    CSVTK_CONCAT(gather(EMMTYPER_MODULE.out.tsv, 'emmtyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = EMMTYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        EMMTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
