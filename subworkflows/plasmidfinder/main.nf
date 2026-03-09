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
 * @output sample_outputs   Per-sample record outputs from PlasmidFinder
 * @output run_outputs    Merged record with consolidated TSV from all samples
 */
nextflow.preview.types = true

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../modules/plasmidfinder/main'
include { CSVTK_CONCAT                          } from '../../modules/csvtk/concat/main'
include { gather                                } from 'plugin/nf-bactopia'

workflow PLASMIDFINDER {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    PLASMIDFINDER_MODULE(assembly)
    CSVTK_CONCAT(gather(PLASMIDFINDER_MODULE.out, 'plasmidfinder', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = PLASMIDFINDER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
