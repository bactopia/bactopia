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
 * @modules csvtk_concat, sccmec as sccmec_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output tsv             Tsv
 * @output merged_tsv      Merged Tsv
 * @output targets         Targets
 * @output target_details  Target Details
 * @output regions         Regions
 * @output regions_details Regions Details
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution logs from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCCMEC {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    SCCMEC_MODULE(fasta)
    CSVTK_CONCAT(gather(SCCMEC_MODULE.out.tsv, 'sccmec'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    targets: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.targets
    target_details: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.target_details
    regions: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.regions
    regions_details: Channel<Tuple<Map, Path>> = SCCMEC_MODULE.out.regions_details

    // Generic aggregate output













    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        SCCMEC_MODULE.out.targets,
        SCCMEC_MODULE.out.target_details,
        SCCMEC_MODULE.out.regions,
        SCCMEC_MODULE.out.regions_details
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SCCMEC_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
