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
 * @output tsv         SISTR comprehensive typing results with serotype and MLST information
 * @output merged_tsv  Combined TSV file containing typing results from all samples
 * @output allele_fasta MLST allele sequences in FASTA format
 * @output allele_json  Detailed MLST allele information in JSON format
 * @output cgmlst_csv   Core genome MLST distance matrix
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SISTR as SISTR_MODULE } from '../../modules/sistr/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow SISTR {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    SISTR_MODULE(assembly)
    CSVTK_CONCAT(gather(SISTR_MODULE.out.tsv, 'sistr'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SISTR_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    allele_fasta: Channel<Tuple<Map, Set<Path>>> = SISTR_MODULE.out.allele_fasta
    allele_json: Channel<Tuple<Map, Set<Path>>> = SISTR_MODULE.out.allele_json
    cgmlst_csv: Channel<Tuple<Map, Set<Path>>> = SISTR_MODULE.out.cgmlst_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        SISTR_MODULE.out.allele_fasta,
        SISTR_MODULE.out.allele_json,
        SISTR_MODULE.out.cgmlst_csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SISTR_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
