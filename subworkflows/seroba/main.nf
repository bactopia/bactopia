/**
 * k-mer based pipeline to identify the serotype of Streptococcus pneumoniae.
 *
 * This subworkflow performs serotyping of Streptococcus pneumoniae from Illumina
 * next-generation sequencing reads using [Seroba](https://github.com/sanger-pathogens/seroba).
 * The tool uses a k-mer based approach to rapidly classify pneumococcal isolates into
 * their respective serotypes based on the capsular polysaccharide synthesis locus.
 *
 * @status stable
 * @keywords Streptococcus, pneumoniae, serotype, k-mer, capsular
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation seroba
 *
 * @modules csvtk_concat, seroba_run
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for S. pneumoniae serotype prediction
 *
 * @output sample_outputs
 * - `tsv`: Serotype prediction results with predicted serotype and confidence in TSV format
 * - `txt`: Detailed information about the predicted serogroup and allele matches
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SEROBA_RUN   } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gather       } from 'plugin/nf-bactopia'

workflow SEROBA {
    take:
    assembly: Channel<Record>

    main:
    SEROBA_RUN(assembly)
    CSVTK_CONCAT(gather(SEROBA_RUN.out, 'seroba', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = SEROBA_RUN.out
    run_outputs = CSVTK_CONCAT.out
}
