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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Bacterial assembly files in FASTA format to be annotated
 *
 * @input proteins
 * Optional protein sequences for homology search. When provided, these trusted protein sequences
 * are used to improve annotation accuracy through homology
 *
 * @input prodigal_tf
 * Optional Prodigal training file. Species-specific training data that improves gene prediction
 * accuracy
 *
 * @output annotations Complete annotation package containing GFF, GBK, FASTA, and other formats
 * @output blastdb        BLAST database created from protein sequences for rapid searching
 * @output faa            Amino acid sequences of predicted proteins in FASTA format
 * @output ffn            Nucleotide sequences of all predicted features (genes, rRNA, tRNA, etc.)
 * @output fna            Input nucleotide assembly (replicated for convenience)
 * @output fsa            Nucleotide sequences of predicted features with pseudogenes
 * @output gbk            GenBank format file suitable for submission to NCBI
 * @output gff            GFF3 format annotation file for use with genome browsers
 * @output tsv            Tab-separated file listing all predicted proteins
 * @output txt            Summary statistics of the annotation results
 * @output sqn            Sequin file for NCBI submission preparation
 * @output tbl            Feature table file for NCBI submission
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'

workflow PROKKA {
    take:
    assembly: Channel<Tuple<Map, Path>>
    proteins: Path
    prodigal_tf: Path?

    main:
    PROKKA_MODULE(assembly, proteins, prodigal_tf)

    emit:
    sample_outputs = PROKKA_MODULE.out
}
