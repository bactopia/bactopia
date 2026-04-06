/**
 * Predict emm types of Streptococcus pyogenes from genome assemblies.
 *
 * This subworkflow uses [emmtyper](https://github.com/MDU-PHL/emmtyper) to predict
 * the emm types of *Streptococcus pyogenes* strains from assembled genomes. It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus pyogenes, emm typing, gas, m protein
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation emmtyper
 *
 * @modules emmtyper, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format to be analyzed for emm genes
 *
 * @input blastdb
 * Optional BLAST database containing emm gene reference sequences for improved typing accuracy
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary of the assigned emm type and cluster
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                 } from 'plugin/nf-bactopia'

workflow EMMTYPER {
    take:
    assembly: Channel<Record>
    blastdb: Value<Path?>

    main:
    EMMTYPER_MODULE(assembly, blastdb)
    CSVTK_CONCAT(gatherCsvtk(EMMTYPER_MODULE.out, 'tsv', [name: 'emmtyper']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = EMMTYPER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
