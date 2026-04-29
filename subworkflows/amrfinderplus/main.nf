/**
 * Find antimicrobial resistance genes and point mutations.
 *
 * This subworkflow uses [AMRFinderPlus](https://github.com/ncbi/amr) to identify acquired antimicrobial
 * resistance genes and some point mutations in protein or assembled nucleotide sequences.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, gene prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,aggregation
 * @citation amrfinderplus
 *
 * @modules csvtk_concat, amrfinderplus_run
 *
 * @input record(meta, genes, proteins, gff)
 * - `meta`: Groovy Record containing sample information
 * - `genes`: Nucleotide sequences of genes in FASTA format
 * - `proteins`: Optional amino acid sequences of proteins in FASTA format (Path?)
 * - `gff`: Optional GFF3 annotation file (Path?)
 *
 * @input db
 * Path to the AMRFinderPlus database directory containing reference data for AMR gene detection.
 *
 * @output sample_outputs
 * - `report`: A tab-delimited report of identified AMR genes and virulence factors
 * - `mutation_report`: Organism-specific point mutations associated with antimicrobial resistance
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with AMRFinder+ results from all samples
 */
nextflow.enable.types = true

include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT      } from '../../modules/csvtk/concat/main'
include { gatherCsvtk       } from 'plugin/nf-bactopia'

workflow AMRFINDERPLUS {
    take:
    fasta: Channel<Record>
    db: Path

    main:
    ch_amrfinderplus_run = AMRFINDERPLUS_RUN(fasta, db)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_amrfinderplus_run, 'report', [name: 'amrfinderplus']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_amrfinderplus_run
    run_outputs = ch_csvtk_concat
}
