/**
 * Predict serogroups of Pseudomonas aeruginosa from assemblies.
 *
 * This subworkflow uses [Pasty](https://github.com/rpetit3/pasty) to perform in silico
 * serogrouping of *Pseudomonas aeruginosa* isolates from assembled genomes. It identifies
 * O-antigen biosynthesis genes to classify isolates into their known serogroups using
 * BLAST-based homology searches.
 *
 * @status stable
 * @keywords pseudomonas aeruginosa, serogroup, typing, o-antigen, prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation pasty
 *
 * @modules pasty, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs   Per-sample record outputs from PASTY_MODULE
 * @output run_outputs    Merged record containing consolidated serogroup predictions from all samples
 */
nextflow.preview.types = true

include { PASTY as PASTY_MODULE } from '../../modules/pasty/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow PASTY {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    PASTY_MODULE(assembly)
    CSVTK_CONCAT(gather(PASTY_MODULE.out, 'pasty', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = PASTY_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
