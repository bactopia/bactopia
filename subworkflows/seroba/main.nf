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
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 * - `tsv`: Serotype prediction results with predicted serotype and confidence in TSV format
 * - `txt`: Detailed information about the predicted serogroup and allele matches
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SEROBA_RUN     } from '../../modules/seroba/run/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { filterWithData } from 'plugin/nf-bactopia'
include { gatherCsvtk         } from 'plugin/nf-bactopia'

workflow SEROBA {
    take:
    reads: Channel<Record>

    main:
    SEROBA_RUN(filterWithData(reads, ['r1', 'r2']))
    CSVTK_CONCAT(gatherCsvtk(SEROBA_RUN.out, 'tsv', [name: 'seroba']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = SEROBA_RUN.out
    run_outputs = CSVTK_CONCAT.out
}
