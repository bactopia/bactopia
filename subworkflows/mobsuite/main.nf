/**
 * Reconstruct and type plasmids from bacterial genome assemblies.
 *
 * This subworkflow uses [MOB-suite](https://github.com/phac-nml/mob-suite) to reconstruct
 * and type plasmids from draft genome assemblies. It separates plasmid from chromosomal
 * sequences, determines plasmid replicon types using the MOB-suite database, and provides
 * comprehensive reports on plasmid content and organization.
 *
 * @status stable
 * @keywords plasmid, reconstruction, typing, mobilome, bacterial genome
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation mobsuite
 *
 * @modules mobsuite_recon, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output txt            Per-sample TSV files containing MOB-suite reconstruction results
 * @output merged_reports Consolidated TSV file containing reconstruction results from all samples
 * @output chromosome      Chromosomal sequences separated from plasmids
 * @output contig_report  Report detailing classification of each contig (chromosome/plasmid)
 * @output plasmids        Reconstructed plasmid sequences
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { flattenPaths   } from 'plugin/nf-bactopia'
include { gather         } from 'plugin/nf-bactopia'

workflow MOBSUITE {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    MOBSUITE_RECON(assembly)
    CSVTK_CONCAT(gather(MOBSUITE_RECON.out.txt, 'mobsuite', 'summary'), 'tsv', 'tsv')

    emit:
    txt: Channel<Tuple<Map, Set<Path>>> = MOBSUITE_RECON.out.txt
    merged_reports: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    chromosome: Channel<Tuple<Map, Set<Path>>> = MOBSUITE_RECON.out.chromosome
    contig_report: Channel<Tuple<Map, Set<Path>>> = MOBSUITE_RECON.out.contig_report
    plasmids: Channel<Tuple<Map, Set<Path>>> = MOBSUITE_RECON.out.plasmids

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
