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
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PROKKA {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    proteins: Path
    prodigal_tf: Path?

    main:
    PROKKA_MODULE(assembly, proteins, prodigal_tf)

    emit:
    // Individual outputs
    annotations: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> = PROKKA_MODULE.out.annotations
    blastdb: Channel<Tuple<Map, Set<Path>>> = PROKKA_MODULE.out.blastdb
    faa: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.faa
    ffn: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.ffn
    fna: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.fna
    fsa: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.fsa
    gbk: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.gbk
    gff: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.gff
    tsv: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.tsv
    txt: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.txt
    sqn: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.sqn
    tbl: Channel<Tuple<Map, Path>> = PROKKA_MODULE.out.tbl

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PROKKA_MODULE.out.blastdb,
        PROKKA_MODULE.out.faa,
        PROKKA_MODULE.out.ffn,
        PROKKA_MODULE.out.fna,
        PROKKA_MODULE.out.fsa,
        PROKKA_MODULE.out.gbk,
        PROKKA_MODULE.out.gff,
        PROKKA_MODULE.out.tsv,
        PROKKA_MODULE.out.txt,
        PROKKA_MODULE.out.sqn,
        PROKKA_MODULE.out.tbl
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PROKKA_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PROKKA_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PROKKA_MODULE.out.versions])
}
