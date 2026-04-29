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
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation busco
 *
 * @modules csvtk_concat, busco
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Genome assemblies to evaluate for completeness. Each record contains metadata
 *   about the sample and a set of assembled contigs in FASTA format.
 *
 * @input busco_lineage
 * BUSCO lineage dataset to use for assessment (e.g., bacteria_odb10).
 *
 * @output sample_outputs
 * - `tsv`: A text summary report of the completeness score (C/S/D/F/M%)
 * - `supplemental`: Directory containing full tables, missing gene lists, and lineage data
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { BUSCO as BUSCO_MODULE } from '../../modules/busco/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk           } from 'plugin/nf-bactopia'

workflow BUSCO {
    take:
    assembly: Channel<Record>
    busco_lineage: String

    main:
    ch_busco = BUSCO_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_busco, 'tsv', [name: "busco-${busco_lineage}"]), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_busco
    run_outputs = ch_csvtk_concat
}
