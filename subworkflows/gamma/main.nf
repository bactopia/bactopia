/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat, gamma as gamma_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @input db
 * Channel containing db data
 *
 * @output gamma        Gamma
 * @output merged_gamma Merged Gamma
 * @output psl          Psl
 * @output fasta        Fasta
 * @output gff          Gff
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution logs from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GAMMA as GAMMA_MODULE } from '../../modules/gamma/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow GAMMA {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    db: Path

    main:
    GAMMA_MODULE(fasta, db)
    CSVTK_CONCAT(gather(GAMMA_MODULE.out.gamma, 'gamma'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gamma: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.gamma
    merged_gamma: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    psl: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.psl
    fasta: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.fasta
    gff: Channel<Tuple<Map, Path>> = GAMMA_MODULE.out.gff

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.gamma,
        CSVTK_CONCAT.out.csv,
        GAMMA_MODULE.out.psl,
        GAMMA_MODULE.out.fasta,
        GAMMA_MODULE.out.gff
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
