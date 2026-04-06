/**
 * Annotate bacterial genomes with functional information.
 *
 * This subworkflow annotates bacterial assemblies using [Prokka](https://github.com/tseemann/prokka).
 * It rapidly calls genes, translates them, and searches them against multiple protein databases
 * to produce comprehensive annotation in various standard formats. Optional protein sequences
 * and Prodigal training files can be provided to improve annotation accuracy.
 *
 * @status stable
 * @keywords bacteria, annotation, genome, prokaryote, functional annotation, genes
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-input
 * @citation prokka
 *
 * @modules prokka as prokka_module
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Bacterial assembly files in FASTA format to be annotated
 *
 * @input proteins
 * Optional protein sequences for homology search to improve annotation accuracy
 *
 * @input prodigal_tf
 * Optional Prodigal training file for improved gene prediction accuracy
 *
 * @output sample_outputs
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbff`: Annotation in GenBank format, containing both sequences and annotations
 * - `fna`: Nucleotide FASTA file of the input contig sequences
 * - `faa`: Protein FASTA file of the translated CDS sequences
 * - `ffn`: Nucleotide FASTA file of all prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
 * - `sqn`: An ASN1 format "Sequin" file for submission to GenBank
 * - `fsa`: Nucleotide FASTA file of the input contig sequences, used by tbl2asn
 * - `tbl`: Feature Table file for NCBI submission
 * - `txt`: Summary statistics relating to the annotated features found
 * - `tsv`: Tab-separated file of all features
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'
include { filterWithData          } from 'plugin/nf-bactopia'

workflow PROKKA {
    take:
    assembly: Channel<Record>
    proteins: Value<Path?>
    prodigal_tf: Value<Path?>

    main:
    PROKKA_MODULE(assembly, proteins, prodigal_tf)

    emit:
    // Downstream inputs
    annotations = filterWithData(PROKKA_MODULE.out, ['fna', 'faa', 'gff'])
    gffs = PROKKA_MODULE.out.map { r -> record(meta: r.meta, gff: r.gff) }
    // Published outputs
    sample_outputs = PROKKA_MODULE.out
    run_outputs = channel.empty()
}
