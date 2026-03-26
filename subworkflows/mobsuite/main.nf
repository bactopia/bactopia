/**
 * Reconstruct and type plasmids from bacterial genome assemblies.
 *
 * This subworkflow uses [MOB-suite](https://github.com/phac-nml/mob-suite) to reconstruct
 * and type plasmids from draft genome assemblies. It separates plasmid from chromosomal
 * sequences, determines plasmid replicon types using the MOB-suite database, and provides
 * comprehensive reports on plasmid content and organization.
 *
 * @status stable
 * @keywords plasmid, reconstruction, typing, mobilome, bacterial genome
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation mobsuite
 *
 * @modules mobsuite_recon, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `chromosome`: Chromosomal sequences separated from plasmid contigs (gzipped FASTA)
 * - `contig_report`: Tab-delimited report assigning each contig to chromosome or plasmid
 * - `txt`: MOB-typer results with replicon type, mobility, and incompatibility group (optional)
 * - `plasmids`: Reconstructed plasmid sequences in gzipped FASTA format (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { gather         } from 'plugin/nf-bactopia'

workflow MOBSUITE {
    take:
    assembly: Channel<Record>

    main:
    MOBSUITE_RECON(assembly)
    CSVTK_CONCAT(gather(MOBSUITE_RECON.out, 'txt', [name: 'mobsuite']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = MOBSUITE_RECON.out
    run_outputs = CSVTK_CONCAT.out
}
