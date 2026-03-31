/**
 * Identify plasmid replicons in bacterial genome assemblies.
 *
 * This subworkflow uses [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder) to
 * identify plasmid replicons in bacterial genome assemblies. It screens assemblies against
 * the PlasmidFinder database to detect known plasmid replicon types and provides detailed
 * results including hit sequences and classification information.
 *
 * @status stable
 * @keywords plasmid, replicon, typing, antimicrobial resistance, mobilome
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation plasmidfinder
 *
 * @modules plasmidfinder, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `json`: PlasmidFinder results in JSON format
 * - `txt`: PlasmidFinder results in text format
 * - `tsv`: Tab-delimited PlasmidFinder results with replicon typing information
 * - `genome_seq`: FASTA sequences of plasmid hits found in the genome (gzipped)
 * - `plasmid_seq`: Reference plasmid sequences matched (gzipped)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../modules/plasmidfinder/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                           } from 'plugin/nf-bactopia'

workflow PLASMIDFINDER {
    take:
    assembly: Channel<Record>

    main:
    PLASMIDFINDER_MODULE(assembly)
    CSVTK_CONCAT(gatherCsvtk(PLASMIDFINDER_MODULE.out, 'tsv', [name: 'plasmidfinder']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = PLASMIDFINDER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
