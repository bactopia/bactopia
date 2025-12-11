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
 * @modules ncbigenomedownload as ncbigenomedownload_module
 *
 * @input accessions
 * Channel containing accessions data
 *
 * @output bactopia_tools Bactopia Tools
 * @output gbk            Gbk
 * @output fna            Fna
 * @output rm             Rm
 * @output features       Features
 * @output gff            Gff
 * @output faa            Faa
 * @output gpff           Gpff
 * @output wgs_gbk        Wgs Gbk
 * @output cds            Cds
 * @output rna            Rna
 * @output rna_fna        Rna Fna
 * @output report         Report
 * @output stats          Stats
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
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
