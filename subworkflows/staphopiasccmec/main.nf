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
 * @output sample_outputs   Per-sample records with SCCmec typing results
 * @output run_outputs    Merged record containing consolidated SCCmec typing from all samples
 */
nextflow.preview.types = true

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE } from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT                              } from '../../modules/csvtk/concat/main'
include { gather                                    } from 'plugin/nf-bactopia'

workflow STAPHOPIASCCMEC {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    STAPHOPIASCCMEC_MODULE(assembly)
    CSVTK_CONCAT(gather(STAPHOPIASCCMEC_MODULE.out, 'staphopiasccmec', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = STAPHOPIASCCMEC_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
