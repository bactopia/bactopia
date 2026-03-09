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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. influenzae serotype prediction
 *
 * @input database_dir
 * Directory containing hicap reference database files (optional)
 *
 * @input model_fp
 * Path to HMM model file for improved detection (optional)
 *
 * @output sample_outputs    Per-sample record outputs from HICAP_MODULE
 * @output run_outputs     Combined hicap results across all samples as a record
 */
nextflow.preview.types = true

include { HICAP as HICAP_MODULE } from '../../modules/hicap/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow HICAP {
    take:
    assembly: Channel<Tuple<Map, Path>>
    database_dir: Path?
    model_fp: Path?

    main:
    HICAP_MODULE(assembly, database_dir, model_fp)
    CSVTK_CONCAT(gather(HICAP_MODULE.out, 'hicap', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = HICAP_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
