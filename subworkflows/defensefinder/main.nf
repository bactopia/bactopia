/**
 * Systematically search for anti-phage defense systems.
 *
 * This subworkflow uses [DefenseFinder](https://github.com/mdmparis/defense-finder) to identify and classify
 * anti-phage defense systems in bacterial genomes. It detects defense genes, HMM hits, and complete
 * defense systems, providing comprehensive analysis of bacterial antiviral mechanisms.
 *
 * @status stable
 * @keywords bacteria, assembly, anti-phage, defense systems, immunity
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation defensefinder
 *
 * @modules defensefinder_run, defensefinder_update, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format for defense system detection
 *
 * @output sample_outputs  Per-sample record outputs from DEFENSEFINDER_RUN
 * @output run_outputs   Combined defense system results across all samples as records
 */
nextflow.preview.types = true

include { DEFENSEFINDER_UPDATE           } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN              } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'
include { gather                         } from 'plugin/nf-bactopia'

workflow DEFENSEFINDER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(assembly, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    GENES_CONCAT(gather(DEFENSEFINDER_RUN.out, 'defensefinder-genes', field: 'genes_tsv'), 'tsv', 'tsv')
    HMMER_CONCAT(gather(DEFENSEFINDER_RUN.out, 'defensefinder-hmmer', field: 'hmmer_tsv'), 'tsv', 'tsv')
    SYSTEMS_CONCAT(gather(DEFENSEFINDER_RUN.out, 'defensefinder-systems', field: 'systems_tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = DEFENSEFINDER_RUN.out
    run_outputs = GENES_CONCAT.out.mix(HMMER_CONCAT.out, SYSTEMS_CONCAT.out)
}
