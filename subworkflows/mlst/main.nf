/**
 * Determine multilocus sequence types (MLST) from bacterial assemblies.
 *
 * This subworkflow uses [mlst](https://github.com/tseemann/mlst) to scan assembled
 * contigs against PubMLST typing schemes and determine sequence types (STs). It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords mlst, sequence typing, pubmlst, bacteria
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation mlst
 *
 * @modules csvtk_concat, mlst
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * PubMLST database to use for MLST typing
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with mlst results from all samples
 */
nextflow.preview.types = true

include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT        } from '../../modules/csvtk/concat/main'
include { gatherCsvtk         } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    assembly: Channel<Record>
    db: Value<Path>

    main:
    ch_mlst = MLST_MODULE(assembly, db)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_mlst, 'tsv', [name: 'mlst']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_mlst
    run_outputs = ch_csvtk_concat
}
