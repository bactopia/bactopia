/**
 * Identify SCCmec elements in Staphylococcus aureus genomes.
 *
 * This subworkflow uses [SCCmec](https://github.com/rpetit3/sccmec) to identify the
 * Staphylococcal Cassette Chromosome mec (SCCmec) element in *Staphylococcus aureus*
 * assemblies. It predicts the type based on the presence of specific *mec* and *ccr*
 * gene complexes, generating detailed BLAST results and typing information.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation sccmec
 *
 * @modules sccmec, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv              Per-sample TSV files with SCCmec typing results
 * @output merged_tsv       Consolidated TSV file containing SCCmec typing from all samples
 * @output targets          Per-sample BLAST results for target sequences
 * @output target_details   Per-sample detailed results for target matches
 * @output regions          Per-sample BLAST results for SCCmec regions
 * @output regions_details  Per-sample detailed results for SCCmec region matches
 * @output results          Aggregated results channel containing all output files
 * @output logs             Aggregated logs channel containing all execution logs
 * @output nf_logs          Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions         Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCCMEC {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    SCCMEC_MODULE(assembly)
    CSVTK_CONCAT(gather(SCCMEC_MODULE.out.tsv, 'sccmec'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Set<Path>>> = SCCMEC_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    targets: Channel<Tuple<Map, Set<Path>>> = SCCMEC_MODULE.out.targets
    target_details: Channel<Tuple<Map, Set<Path>>> = SCCMEC_MODULE.out.target_details
    regions: Channel<Tuple<Map, Set<Path>>> = SCCMEC_MODULE.out.regions
    regions_details: Channel<Tuple<Map, Set<Path>>> = SCCMEC_MODULE.out.regions_details

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
