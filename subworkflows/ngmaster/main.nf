/**
 * Perform multi-antigen sequence typing of Neisseria gonorrhoeae from genome assemblies.
 *
 * This subworkflow uses [ngmaster](https://github.com/MDU-PHL/ngmaster) to perform
 * in silico multi-antigen sequence typing (NG-MAST) for *Neisseria gonorrhoeae*
 * strains from assembled genomes. It processes each sample individually and aggregates
 * the results into a single consolidated report.
 *
 * @status stable
 * @keywords neisseria gonorrhoeae, ng-mast, typing, gonococcal, antigen
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation ngmaster
 *
 * @modules csvtk_concat, ngmaster
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited NG-MASTER results with porB and tbpB alleles and sequence type
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { NGMASTER as NGMASTER_MODULE } from '../../modules/ngmaster/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                 } from 'plugin/nf-bactopia'

workflow NGMASTER {
    take:
    assembly: Channel<Record>

    main:
    ch_ngmaster = NGMASTER_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_ngmaster, 'tsv', [name: 'ngmaster']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_ngmaster
    run_outputs = ch_csvtk_concat
}
