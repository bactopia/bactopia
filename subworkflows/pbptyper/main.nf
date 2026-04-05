/**
 * Predict penicillin binding protein (PBP) types of Streptococcus pneumoniae from genome assemblies.
 *
 * This subworkflow uses [pbptyper](https://github.com/rpetit3/pbptyper) to predict
 * the penicillin binding protein (PBP) types and predict antimicrobial susceptibility
 * of *Streptococcus pneumoniae* strains from assembled genomes. It processes each sample
 * individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, pbp typing, penicillin, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation pbptyper
 *
 * @modules pbptyper, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary file with the predicted PBP type for each gene
 * - `blast`: A tab-delimited file of the raw TBLASTN hits used for gene identification
 * - `details`: Detailed PBP typing results for each gene analyzed
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                 } from 'plugin/nf-bactopia'

workflow PBPTYPER {
    take:
    assembly: Channel<Record>

    main:
    PBPTYPER_MODULE(assembly)
    CSVTK_CONCAT(gatherCsvtk(PBPTYPER_MODULE.out, 'tsv', [name: 'pbptyper']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = PBPTYPER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
