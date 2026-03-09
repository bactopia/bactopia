/**
 * Determine multilocus sequence types (MLST) from bacterial assemblies.
 *
 * This subworkflow uses [mlst](https://github.com/tseemann/mlst) to scan assembled
 * contigs against PubMLST typing schemes and determine sequence types (STs). It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords mlst, sequence typing, pubmlst, bacteria
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation mlst
 *
 * @modules csvtk_concat, mlst
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * PubMLST database to use for MLST typing
 *
 * @output sample_outputs   Per-sample records with MLST typing results
 * @output run_outputs    Merged record containing consolidated MLST typing from all samples
 */
nextflow.preview.types = true

include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT        } from '../../modules/csvtk/concat/main'
include { gather              } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    assembly: Channel<Tuple<Map, Path>>
    db: Path

    main:
    MLST_MODULE(assembly, db)
    CSVTK_CONCAT(gather(MLST_MODULE.out, 'mlst', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = MLST_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
