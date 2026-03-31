/**
 * In silico serotyping of the Haemophilus influenzae capsule locus.
 *
 * This subworkflow performs serotyping of Haemophilus influenzae by analyzing
 * the capsule locus using [hicap](https://github.com/scwatts/hicap). The tool
 * identifies capsule type and produces detailed annotations of the cap locus,
 * including graphical representation in SVG format and annotation in GenBank format.
 *
 * @status stable
 * @keywords Haemophilus, influenzae, serotype, capsule, cap locus
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-input, aggregation
 * @citation hicap
 *
 * @modules csvtk_concat, hicap
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. influenzae serotype prediction
 *
 * @input database_dir
 * Directory containing hicap reference database files (optional)
 *
 * @input model_fp
 * Path to HMM model file for improved detection (optional)
 *
 * @output sample_outputs
 * - `gbff`: GenBank file containing the annotated capsule locus region (optional)
 * - `svg`: SVG visualization of the capsule locus gene arrangement (optional)
 * - `tsv`: Tab-delimited summary of the predicted serotype and locus coverage
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk           } from 'plugin/nf-bactopia'

workflow HICAP {
    take:
    assembly: Channel<Record>
    database_dir: Path?
    model_fp: Path?

    main:
    HICAP_MODULE(assembly, database_dir, model_fp)
    CSVTK_CONCAT(gatherCsvtk(HICAP_MODULE.out, 'tsv', [name: 'hicap']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = HICAP_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
