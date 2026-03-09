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
 * @output sample_outputs   Per-sample records with SCCmec typing results
 * @output run_outputs    Merged record containing consolidated SCCmec typing from all samples
 */
nextflow.preview.types = true

include { SCCMEC as SCCMEC_MODULE } from '../../modules/sccmec/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCCMEC {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    SCCMEC_MODULE(assembly)
    CSVTK_CONCAT(gather(SCCMEC_MODULE.out, 'sccmec', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = SCCMEC_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
