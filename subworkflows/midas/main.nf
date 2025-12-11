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
 * @modules csvtk_concat, midas_species
 *
 * @input reads
 * Channel containing reads data
 *
 * @input database
 * Channel containing database data
 *
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output abundances Abundances
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MIDAS_SPECIES } from '../../modules/midas/species/main'
include { CSVTK_CONCAT  } from '../../modules/csvtk/concat/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow MIDAS {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    MIDAS_SPECIES(reads, database)
    CSVTK_CONCAT(gather(MIDAS_SPECIES.out.tsv, 'midas'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MIDAS_SPECIES.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    abundances: Channel<Tuple<Map, Path>> = MIDAS_SPECIES.out.abundances

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.tsv,
        MIDAS_SPECIES.out.abundances,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MIDAS_SPECIES.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
