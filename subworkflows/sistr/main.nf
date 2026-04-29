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
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembly files in FASTA format for Salmonella typing
 *
 * @output sample_outputs
 * - `tsv`: SISTR prediction results in TSV format
 * - `allele_fasta`: Novel alleles in FASTA format
 * - `allele_json`: Alleles in JSON format
 * - `cgmlst_csv`: cgMLST profile in CSV format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { SISTR as SISTR_MODULE } from '../../modules/sistr/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk           } from 'plugin/nf-bactopia'

workflow SISTR {
    take:
    assembly: Channel<Record>

    main:
    ch_sistr = SISTR_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_sistr, 'tsv', [name: 'sistr']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_sistr
    run_outputs = ch_csvtk_concat
}
