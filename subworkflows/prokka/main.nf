/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules prokka as prokka_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @input proteins
 * Channel containing proteins data
 *
 * @input prodigal_tf
 * Channel containing prodigal_tf data
 *
 * @output annotations Annotations
 * @output blastdb     Blastdb
 * @output faa         Faa
 * @output ffn         Ffn
 * @output fna         Fna
 * @output fsa         Fsa
 * @output gbk         Gbk
 * @output gff         Gff
 * @output tsv         Tsv
 * @output txt         Txt
 * @output sqn         Sqn
 * @output tbl         Tbl
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution logs from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PROKKA as PROKKA_MODULE } from '../../modules/prokka/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PROKKA {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>
    proteins: Path
    prodigal_tf: Path?

    main:
    PROKKA_MODULE(fasta, proteins, prodigal_tf)

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
