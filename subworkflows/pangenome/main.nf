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
 * @subworkflows pirate, roary, panaroo, snpdists
 *
 * @input gff
 * Channel containing gff data
 *
 * @input use_pirate
 * Channel containing use_pirate data
 *
 * @input use_roary
 * Channel containing use_roary data
 *
 * @output aln      Aln
 * @output csv      Csv
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PIRATE       } from '../pirate/main'
include { ROARY        } from '../roary/main'
include { PANAROO      } from '../panaroo/main'
include { SNPDISTS     } from '../snpdists/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow PANGENOME {
    take:
    gff        : Channel<Tuple<Map, Set<Path>>>
    use_pirate : Boolean
    use_roary  : Boolean

    main:

    // Initialize channels
    ch_aln = channel.empty() as Channel<Tuple<Map, Path>>
    ch_csv = channel.empty() as Channel<Tuple<Map, Path>>
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    // Choose pangenome tool based on params
    if (use_pirate) {
        PIRATE(gff)
        ch_aln = PIRATE.out.aln
        ch_csv = PIRATE.out.csv
        ch_results = PIRATE.out.results
        ch_logs = PIRATE.out.logs
        ch_nf_logs = PIRATE.out.nf_logs
        ch_versions = PIRATE.out.versions
    } else if (use_roary) {
        ROARY(gff)
        ch_aln = ROARY.out.aln
        ch_csv = ROARY.out.csv
        ch_results = ROARY.out.results
        ch_logs = ROARY.out.logs
        ch_nf_logs = ROARY.out.nf_logs
        ch_versions = ROARY.out.versions
    } else {
        PANAROO(gff)
        ch_aln = PANAROO.out.filtered_aln
        ch_csv = PANAROO.out.csv
        ch_results = PANAROO.out.results
        ch_logs = PANAROO.out.logs
        ch_nf_logs = PANAROO.out.nf_logs
        ch_versions = PANAROO.out.versions
    }

    // Per-sample SNP distances
    ch_unmasked_aln = ch_aln.map({ _meta, aln -> 
        tuple([name: "core-genome.distance", process_name: "snpdists"], aln)
    })
    SNPDISTS(ch_unmasked_aln)
    ch_results = ch_results.mix(SNPDISTS.out.results)
    ch_logs = ch_logs.mix(SNPDISTS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNPDISTS.out.nf_logs)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Path>> = ch_aln
    csv: Channel<Tuple<Map, Path>> = ch_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([ch_results])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ch_versions])
}
