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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation pbptyper
 *
 * @modules pbptyper, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs   Per-sample record outputs from PBPTYPER_MODULE
 * @output run_outputs    Merged record with consolidated TSV from all samples
 */
nextflow.preview.types = true

include { PBPTYPER as PBPTYPER_MODULE } from '../../modules/pbptyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { gather                      } from 'plugin/nf-bactopia'

workflow PBPTYPER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    PBPTYPER_MODULE(assembly)
    CSVTK_CONCAT(gather(PBPTYPER_MODULE.out, 'pbptyper', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = PBPTYPER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
