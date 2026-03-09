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
 * @output sample_outputs  Per-sample records containing meta, summary, results, logs, nf_logs, versions
 * @output run_outputs   Cross-sample aggregation record
 */
nextflow.preview.types = true

include { AGRVATE as AGRVATE_MODULE } from '../../modules/agrvate/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gather                    } from 'plugin/nf-bactopia'

workflow AGRVATE {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    AGRVATE_MODULE(assembly)
    CSVTK_CONCAT(gather(AGRVATE_MODULE.out, 'agrvate', field: 'summary'), 'tsv', 'tsv')

    emit:
    // Per-sample records (contains meta, summary, results, logs, nf_logs, versions)
    sample_outputs = AGRVATE_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
