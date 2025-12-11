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
 * @modules sketcher as sketcher_module
 *
 * @input reads
 * Channel containing reads data
 *
 * @input mash_db
 * Channel containing mash_db data
 *
 * @input sourmash_db
 * Channel containing sourmash_db data
 *
 * @output sig      Sig
 * @output msh      Msh
 * @output mash     Mash
 * @output sourmash Sourmash
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SKETCHER as SKETCHER_MODULE } from '../../../modules/bactopia/sketcher/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SKETCHER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    mash_db: Path
    sourmash_db: Path

    main:
    SKETCHER_MODULE(reads, mash_db, sourmash_db)

    emit:
    // Individual outputs
    sig: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sig
    msh: Channel<Tuple<Map, Set<Path>>> = SKETCHER_MODULE.out.msh
    mash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.mash
    sourmash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sourmash

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SKETCHER_MODULE.out.sig,
        SKETCHER_MODULE.out.msh,
        SKETCHER_MODULE.out.mash,
        SKETCHER_MODULE.out.sourmash
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.versions])
}
