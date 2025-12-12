/**
 * Download bacterial genomes from NCBI's RefSeq database.
 *
 * This subworkflow downloads complete and draft bacterial genomes using the
 * [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) tool. It fetches
 * genome assemblies in various formats including GenBank, GFF, and FASTA files
 * along with associated annotation files and statistics.
 *
 * @status stable
 * @keywords download, ncbi, refseq, genome, assembly, database
 * @tags complexity:moderate input-type:single output-type:multiple features:resource-download database-dependent
 * @citation ncbigenomedownload
 *
 * @modules ncbigenomedownload as ncbigenomedownload_module
 *
 * @input accessions
 * A file containing NCBI accession numbers, one per line. If empty, will download all genomes matching the specified criteria.
 *
 * @output bactopia_tools A channel of downloaded files formatted for Bactopia Tools workflows
 * @output gbk              GenBank format files containing complete genome annotations
 * @output fna              FASTA format files containing genomic sequences
 * @output rm               README.md files with metadata about downloaded genomes
 * @output features         Feature table files in text format
 * @output gff              GFF3 format files with genome feature annotations
 * @output faa              Protein FASTA files of translated CDS sequences
 * @output gpff             Genomic protein FASTA files
 * @output wgs_gbk          Whole Genome Shotgun (WGS) project GenBank files
 * @output cds              Nucleotide FASTA files of CDS sequences
 * @output rna              RNA FASTA files (rRNA, tRNA, other RNA)
 * @output rna_fna          RNA sequences in FASTA format
 * @output report           Summary report of downloaded genomes
 * @output stats            Statistics files with assembly metrics
 * @output results          Aggregated results channel containing all output files
 * @output logs             Aggregated logs channel containing all execution logs
 * @output nf_logs          Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions         Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../modules/ncbigenomedownload/main'
include { flattenPaths                                    } from 'plugin/nf-bactopia'
include { gather                                          } from 'plugin/nf-bactopia'

workflow NCBIGENOMEDOWNLOAD {

    take:
    accessions: Path?

    main:
    NCBIGENOMEDOWNLOAD_MODULE(accessions)
    ch_to_bactopia_tools = NCBIGENOMEDOWNLOAD_MODULE.out.all.map { path -> [[id: file(path).getSimpleName()], file(path)] }

    emit:
    // Individual outputs
    bactopia_tools: Channel<Tuple<Map, Path>> = ch_to_bactopia_tools
    gbk: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gbk
    fna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.fna
    rm: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rm
    features: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.features
    gff: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gff
    faa: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.faa
    gpff: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.gpff
    wgs_gbk: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk
    cds: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.cds
    rna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rna
    rna_fna: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna
    report: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.report
    stats: Channel<Tuple<Map, Path>> = NCBIGENOMEDOWNLOAD_MODULE.out.stats

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        NCBIGENOMEDOWNLOAD_MODULE.out.gbk,
        NCBIGENOMEDOWNLOAD_MODULE.out.fna,
        NCBIGENOMEDOWNLOAD_MODULE.out.rm,
        NCBIGENOMEDOWNLOAD_MODULE.out.features,
        NCBIGENOMEDOWNLOAD_MODULE.out.gff,
        NCBIGENOMEDOWNLOAD_MODULE.out.faa,
        NCBIGENOMEDOWNLOAD_MODULE.out.gpff,
        NCBIGENOMEDOWNLOAD_MODULE.out.wgs_gbk,
        NCBIGENOMEDOWNLOAD_MODULE.out.cds,
        NCBIGENOMEDOWNLOAD_MODULE.out.rna,
        NCBIGENOMEDOWNLOAD_MODULE.out.rna_fna,
        NCBIGENOMEDOWNLOAD_MODULE.out.report,
        NCBIGENOMEDOWNLOAD_MODULE.out.stats
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([NCBIGENOMEDOWNLOAD_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([NCBIGENOMEDOWNLOAD_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([NCBIGENOMEDOWNLOAD_MODULE.out.versions])
}
