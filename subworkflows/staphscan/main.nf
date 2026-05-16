/**
 * Genome-based surveillance analysis of Staphylococcus aureus.
 *
 * This subworkflow uses [StaphSCAN](https://github.com/riccabolla/StaphSCAN) to perform
 * genome-based surveillance of *Staphylococcus aureus*, integrating species identification,
 * MLST, *spa* typing, SCCmec typing, capsular typing, and detection of virulence, biofilm,
 * and antimicrobial resistance genes. It processes each sample individually and aggregates
 * the results into a single consolidated report.
 *
 * @status stable
 * @keywords staphylococcus aureus, surveillance, mlst, spa typing, sccmec, amr, virulence
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation staphscan
 *
 * @modules csvtk_concat, staphscan
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * Custom MLST database directory
 *
 * @output sample_outputs
 * - `tsv`: Per-sample surveillance summary with MLST, spa type, SCCmec, capsule, AGR, resistance, biofilm, and virulence results
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with staphscan results from all samples
 */
nextflow.enable.types = true

include { STAPHSCAN as STAPHSCAN_MODULE } from '../../modules/staphscan/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow STAPHSCAN {
    take:
    fna: Channel<Record>
    db: Path?

    main:
    ch_staphscan = STAPHSCAN_MODULE(fna, db)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_staphscan, 'tsv', [name: 'staphscan']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_staphscan
    run_outputs = ch_csvtk_concat
}
