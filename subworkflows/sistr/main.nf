/**
 * Salmonella In Silico Typing Resource command-line tool.
 *
 * This subworkflow performs comprehensive typing of Salmonella genomes using
 * [SISTR](https://github.com/phac-nml/sistr_cmd), which predicts serotype,
 * determines subspecies, performs MLST typing, and calculates core genome
 * MLST distances. The tool provides a one-stop solution for Salmonella
 * classification and epidemiological typing.
 *
 * @status stable
 * @keywords Salmonella, serotype, MLST, cgMLST, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation sistr
 *
 * @modules csvtk_concat, sistr
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for Salmonella typing
 *
 * @output sample_outputs  Record outputs from SISTR module per sample
 * @output run_outputs   Record with merged CSV and associated logs/versions
 */
nextflow.preview.types = true

include { SISTR as SISTR_MODULE } from '../../modules/sistr/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow SISTR {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    SISTR_MODULE(assembly)
    CSVTK_CONCAT(gather(SISTR_MODULE.out, 'sistr', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = SISTR_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
