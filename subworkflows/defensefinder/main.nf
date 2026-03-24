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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format for defense system detection
 *
 * @output sample_outputs
 * - `genes_tsv`: Tab-delimited list of detected defense genes
 * - `hmmer_tsv`: Tab-delimited list of HMMER hits used for detection
 * - `systems_tsv`: Tab-delimited summary of detected defense systems
 * - `proteins`: Protein sequences of the detected defense genes
 * - `proteins_index`: Index file for the protein sequences
 * - `macsydata_raw`: Compressed tarball of raw MacSyFinder data (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
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
    assembly: Channel<Record>

    main:
    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(assembly, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    GENES_CONCAT(gather(DEFENSEFINDER_RUN.out, 'genes_tsv', [name: 'defensefinder-genes']), 'tsv', 'tsv')
    HMMER_CONCAT(gather(DEFENSEFINDER_RUN.out, 'hmmer_tsv', [name: 'defensefinder-hmmer']), 'tsv', 'tsv')
    SYSTEMS_CONCAT(gather(DEFENSEFINDER_RUN.out, 'systems_tsv', [name: 'defensefinder-systems']), 'tsv', 'tsv')

    emit:
    sample_outputs = DEFENSEFINDER_RUN.out
    run_outputs = GENES_CONCAT.out.mix(HMMER_CONCAT.out, SYSTEMS_CONCAT.out)
}
