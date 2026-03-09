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
 * @output sample_outputs    Per-sample record outputs from BUSCO_MODULE
 * @output run_outputs     Combined BUSCO results across all samples as a record
 */
nextflow.preview.types = true

include { BUSCO as BUSCO_MODULE } from '../../modules/busco/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow BUSCO {
    take:
    assembly: Channel<Tuple<Map, Path>>
    busco_lineage: String

    main:
    BUSCO_MODULE(assembly)
    CSVTK_CONCAT(gather(BUSCO_MODULE.out, "busco-${busco_lineage}", field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = BUSCO_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
