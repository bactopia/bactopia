/**
 * Systematically search for anti-phage defense systems.
 *
 * This subworkflow uses [DefenseFinder](https://github.com/mdmparis/defense-finder) to identify and classify
 * anti-phage defense systems in bacterial genomes. It detects defense genes, HMM hits, and complete
 * defense systems, providing comprehensive analysis of bacterial antiviral mechanisms.
 *
 * @status stable
 * @keywords bacteria, assembly, anti-phage, defense systems, immunity
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation defensefinder
 *
 * @modules defensefinder_run, defensefinder_update, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
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
// bactopia-lint: ignore S015
nextflow.preview.types = true

include { DEFENSEFINDER_UPDATE           } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN              } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                    } from 'plugin/nf-bactopia'

workflow DEFENSEFINDER {
    take:
    assembly: Channel<Record>

    main:
    ch_defensefinder_update = DEFENSEFINDER_UPDATE()
    ch_defensefinder_run = DEFENSEFINDER_RUN(assembly, ch_defensefinder_update.map { r -> r.db })

    // Merge results
    ch_genes_concat = GENES_CONCAT(gatherCsvtk(ch_defensefinder_run, 'genes_tsv', [name: 'defensefinder-genes']), 'tsv', 'tsv')
    ch_hmmer_concat = HMMER_CONCAT(gatherCsvtk(ch_defensefinder_run, 'hmmer_tsv', [name: 'defensefinder-hmmer']), 'tsv', 'tsv')
    ch_systems_concat = SYSTEMS_CONCAT(gatherCsvtk(ch_defensefinder_run, 'systems_tsv', [name: 'defensefinder-systems']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_defensefinder_run
    run_outputs = ch_genes_concat.mix(ch_hmmer_concat, ch_systems_concat)
}
