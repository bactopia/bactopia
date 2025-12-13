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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv          Per-sample TSV files containing plasmid replicon detection results
 * @output merged_tsv   Consolidated TSV file containing plasmid replicon results from all samples
 * @output json         JSON files with detailed plasmid replicon detection information
 * @output txt          Summary text files with plasmid replicon findings
 * @output genome_seq   FASTA files containing genome sequences with plasmid hits
 * @output plasmid_seq  FASTA files containing identified plasmid sequences
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../modules/plasmidfinder/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { flattenPaths                          } from 'plugin/nf-bactopia'
include { gather                                } from 'plugin/nf-bactopia'

workflow PLASMIDFINDER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    PLASMIDFINDER_MODULE(assembly)
    CSVTK_CONCAT(gather(PLASMIDFINDER_MODULE.out.tsv, 'plasmidfinder'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = PLASMIDFINDER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    json: Channel<Tuple<Map, Set<Path>>> = PLASMIDFINDER_MODULE.out.json
    txt: Channel<Tuple<Map, Set<Path>>> = PLASMIDFINDER_MODULE.out.txt
    genome_seq: Channel<Tuple<Map, Set<Path>>> = PLASMIDFINDER_MODULE.out.genome_seq
    plasmid_seq: Channel<Tuple<Map, Set<Path>>> = PLASMIDFINDER_MODULE.out.plasmid_seq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PLASMIDFINDER_MODULE.out.json,
        PLASMIDFINDER_MODULE.out.txt,
        PLASMIDFINDER_MODULE.out.genome_seq,
        PLASMIDFINDER_MODULE.out.plasmid_seq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PLASMIDFINDER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
