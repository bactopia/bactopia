/**
 * Identify SCCmec elements in Staphylococcus aureus genomes using Staphopia method.
 *
 * This subworkflow uses [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec) to
 * identify Staphylococcal Cassette Chromosome mec (SCCmec) elements in *Staphylococcus aureus*
 * assemblies. This is the standalone version of the SCCmec typing method developed for
 * the Staphopia project, which predicts SCCmec types based on the presence of specific
 * *mec* and *ccr* gene complexes.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation staphopiasccmec
 *
 * @modules staphopiasccmec, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing SCCmec typing results
 * @output merged_tsv  Consolidated TSV file containing SCCmec typing results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE } from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT                              } from '../../modules/csvtk/concat/main'
include { flattenPaths                              } from 'plugin/nf-bactopia'
include { gather                                    } from 'plugin/nf-bactopia'

workflow STAPHOPIASCCMEC {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    STAPHOPIASCCMEC_MODULE(assembly)
    CSVTK_CONCAT(gather(STAPHOPIASCCMEC_MODULE.out.tsv, 'staphopiasccmec'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
