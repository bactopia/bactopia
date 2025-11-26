//
// mobsuite - Reconstruct and annotate plasmids in bacterial assemblies
//
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow MOBSUITE {
    take:
    fasta: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ assemblies ] ]

    main:
    MOBSUITE_RECON(fasta)
    CSVTK_CONCAT(gather(MOBSUITE_RECON.out.txt, 'mobsuite', 'summary'), 'tsv', 'tsv')

    emit:
    txt: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.txt
    merged_reports: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    chromosome: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.chromosome
    contig_report: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.contig_report
    plasmids: Channel<Tuple<Map, Path>> = MOBSUITE_RECON.out.plasmids

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.txt,
        CSVTK_CONCAT.out.csv,
        MOBSUITE_RECON.out.chromosome,
        MOBSUITE_RECON.out.contig_report,
        MOBSUITE_RECON.out.plasmids
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MOBSUITE_RECON.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
