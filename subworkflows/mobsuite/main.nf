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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs  Per-sample records containing meta, chromosome, contig_report, txt, plasmids, results, logs, nf_logs, and versions
 * @output run_outputs   Merged record containing meta, csv, results, logs, nf_logs, and versions
 */
nextflow.preview.types = true

include { MOBSUITE_RECON } from '../../modules/mobsuite/recon/main'
include { CSVTK_CONCAT   } from '../../modules/csvtk/concat/main'
include { gather         } from 'plugin/nf-bactopia'

workflow MOBSUITE {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    MOBSUITE_RECON(assembly)
    CSVTK_CONCAT(gather(MOBSUITE_RECON.out, 'mobsuite', field: 'txt'), 'tsv', 'tsv')

    emit:
    sample_outputs = MOBSUITE_RECON.out
    run_outputs = CSVTK_CONCAT.out
}
