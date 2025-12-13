/**
 * Identify Staphylococcus aureus agr locus type and operon variants.
 *
 * This subworkflow uses [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) to rapidly identify the
 * accessory gene regulator (agr) locus type and detect agr operon variants in Staphylococcus aureus.
 * The agr system is a key quorum-sensing regulator of virulence in S. aureus.
 *
 * @status stable
 * @keywords staphylococcus aureus, assembly, agr locus, virulence, quorum sensing
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation agrvate
 *
 * @modules agrvate, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format for agr locus detection
 *
 * @output tsv          Agr locus typing results in TSV format for each sample
 * @output supplemental Additional detailed results including variant analysis
 * @output merged_tsv   Combined agr typing results from all samples in a single TSV file
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { AGRVATE as AGRVATE_MODULE } from '../../modules/agrvate/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow AGRVATE {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    AGRVATE_MODULE(assembly)
    CSVTK_CONCAT(gather(AGRVATE_MODULE.out.summary, 'agrvate'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Set<Path>>> = AGRVATE_MODULE.out.summary
    supplemental: Channel<Tuple<Map, Set<Path>>> = AGRVATE_MODULE.out.supplemental
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate output
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.summary,
        AGRVATE_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AGRVATE_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
