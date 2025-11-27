//
// panaroo - Pipeline for pangenome investigations
//
nextflow.preview.types = true

include { PANAROO_RUN  } from '../../modules/panaroo/run/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow PANAROO {
    take:
    gff: Channel<Tuple<Map, Set<Path>>>

    main:
    PANAROO_RUN(gff)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.csv
    aln: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.aln
    filtered_aln: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.filtered_aln
    panaroo_csv: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.panaroo_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PANAROO_RUN.out.supplemental,
        PANAROO_RUN.out.csv,
        PANAROO_RUN.out.aln,
        PANAROO_RUN.out.filtered_aln,
        PANAROO_RUN.out.panaroo_csv
    ])
    logs: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.versions
}
