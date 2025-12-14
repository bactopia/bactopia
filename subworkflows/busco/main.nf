/**
 * Assess genome assembly completeness using BUSCO.
 *
 * This subworkflow evaluates genome assembly completeness by searching for
 * single-copy orthologs against the [BUSCO](https://busco.ezlab.org/) database.
 * It generates comprehensive completeness reports including missing, duplicated,
 * fragmented, and complete single-copy orthologs. The workflow includes individual
 * sample assessments and a merged summary report across all samples.
 *
 * @status stable
 * @keywords assembly, completeness, quality, assessment, orthologs, evaluation
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation busco
 *
 * @modules csvtk_concat, busco
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Genome assemblies to evaluate for completeness. Each tuple contains metadata
 *   about the sample and a set of assembled contigs in FASTA format.
 *
 * @input busco_lineage
 * BUSCO lineage dataset to use for assessment (e.g., bacteria_odb10).
 *
 * @output tsv          Per-sample BUSCO completeness assessment results in TSV format
 * @output merged_tsv   Combined BUSCO results summary across all samples
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BUSCO as BUSCO_MODULE } from '../../modules/busco/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow BUSCO {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    busco_lineage: String

    main:
    BUSCO_MODULE(assembly)
    CSVTK_CONCAT(gather(BUSCO_MODULE.out.tsv, "busco-${busco_lineage}"), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = BUSCO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.tsv,
        BUSCO_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BUSCO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
