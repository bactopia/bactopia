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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation ngmaster
 *
 * @modules csvtk_concat, ngmaster
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs  Per-sample records containing NG-MAST typing results, logs, and versions
 * @output run_outputs   Merged record containing consolidated NG-MAST typing from all samples
 */
nextflow.preview.types = true

include { NGMASTER as NGMASTER_MODULE } from '../../modules/ngmaster/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gather                      } from 'plugin/nf-bactopia'

workflow NGMASTER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    NGMASTER_MODULE(assembly)
    CSVTK_CONCAT(gather(NGMASTER_MODULE.out, 'ngmaster', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = NGMASTER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
