/**
 * Predict serotypes of Shigella and EIEC from assemblies.
 *
 * This subworkflow uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to predict
 * serotypes of *Shigella* and Enteroinvasive *E. coli* (EIEC) from assembled genomes.
 * It uses a cluster-informed approach to identify specific serotype markers and classify
 * isolates based on their antigenic profiles.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, typing, cluster analysis
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation shigeifinder
 *
 * @modules shigeifinder, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs   Per-sample records containing ShigEiFinder serotype predictions
 * @output run_outputs    Merged record containing consolidated serotype predictions from all samples
 */
nextflow.preview.types = true

include { SHIGEIFINDER as SHIGEIFINDER_MODULE } from '../../modules/shigeifinder/main'
include { CSVTK_CONCAT                        } from '../../modules/csvtk/concat/main'
include { gather                              } from 'plugin/nf-bactopia'

workflow SHIGEIFINDER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    SHIGEIFINDER_MODULE(assembly)
    CSVTK_CONCAT(gather(SHIGEIFINDER_MODULE.out, 'shigeifinder', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = SHIGEIFINDER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
