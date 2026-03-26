/**
 * Prediction of prophages from bacterial genomes.
 *
 * This subworkflow identifies prophages in bacterial genomes using [PhiSpy](https://github.com/linsalrob/PhiSpy),
 * which combines similarity-based and composition-based strategies for accurate detection.
 * The tool identifies integrated phage sequences, extracts bacterial and phage regions,
 * and provides comprehensive annotation including GFF format for downstream analysis.
 *
 * @status stable
 * @keywords prophage, phage, bacterial, genome, mobile genetic elements
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation phispy
 *
 * @modules csvtk_concat, phispy
 *
 * @input gbk
 * Annotated bacterial genomes in GenBank format for prophage prediction
 *
 * @output sample_outputs
 * - `tsv`: Coordinates (start/end) of each predicted prophage region in the genome
 * - `supplemental`: Directory containing detailed prophage information, sequences, and annotations
 *
 * @output run_outputs
 * - `csv`: Merged prophage prediction results from all samples
 */
nextflow.preview.types = true

include { PHISPY as PHISPY_MODULE } from '../../modules/phispy/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gather                  } from 'plugin/nf-bactopia'

workflow PHISPY {
    take:
    gbk: Channel<Record>

    main:
    PHISPY_MODULE(gbk)
    CSVTK_CONCAT(gather(PHISPY_MODULE.out, 'tsv', [name: 'phispy']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = PHISPY_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
