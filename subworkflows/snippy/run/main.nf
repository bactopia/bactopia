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
 * @modules snippy_run
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Sequencing reads for variant calling
 *
 * @input reference
 * Channel containing reference data
 *
 * @output aligned_fa               Aligned Fa
 * @output vcf                      Vcf
 * @output aligned_fa_error         Aligned Fa Error
 * @output vcf_error                Vcf Error
 * @output error                    Error
 * @output annotated_vcf            Annotated Vcf
 * @output bam                      Bam
 * @output bai                      Bai
 * @output bed                      Bed
 * @output consensus_fa             Consensus Fa
 * @output consensus_subs_fa        Consensus Subs Fa
 * @output consensus_subs_masked_fa Consensus Subs Masked Fa
 * @output coverage                 Coverage
 * @output csv                      Csv
 * @output filt_vcf                 Filt Vcf
 * @output gff                      Gff
 * @output html                     Html
 * @output raw_vcf                  Raw Vcf
 * @output subs_vcf                 Subs Vcf
 * @output tab                      Tab
 * @output txt                      Txt
 * @output results                  Aggregated results channel containing all output files
 * @output logs                     Aggregated logs channel containing all execution logs
 * @output nf_logs                  Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions                 Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SNIPPY_RUN   } from '../../../modules/snippy/run/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow SNIPPY {
    take:
    reads : Channel<Tuple<Map, Set<Path>>>
    reference : Path

    main:
    SNIPPY_RUN(reads, reference)

    emit:
    aligned_fa: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.aligned_fa
    vcf: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.vcf
    aligned_fa_error: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.aligned_fa_error
    vcf_error: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.vcf_error
    error: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.error
    annotated_vcf: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.annotated_vcf
    bam: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.bam
    bai: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.bai
    bed: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.bed
    consensus_fa: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.consensus_fa
    consensus_subs_fa: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.consensus_subs_fa
    consensus_subs_masked_fa: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.consensus_subs_masked_fa
    coverage: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.coverage
    csv: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.csv
    filt_vcf: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.filt_vcf
    gff: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.gff
    html: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.html
    raw_vcf: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.raw_vcf
    subs_vcf: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.subs_vcf
    tab: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.tab
    txt: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.txt

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SNIPPY_RUN.out.aligned_fa,
        SNIPPY_RUN.out.vcf,
        SNIPPY_RUN.out.annotated_vcf,
        SNIPPY_RUN.out.bam,
        SNIPPY_RUN.out.bai,
        SNIPPY_RUN.out.bed,
        SNIPPY_RUN.out.consensus_fa,
        SNIPPY_RUN.out.consensus_subs_fa,
        SNIPPY_RUN.out.consensus_subs_masked_fa,
        SNIPPY_RUN.out.coverage,
        SNIPPY_RUN.out.csv,
        SNIPPY_RUN.out.filt_vcf,
        SNIPPY_RUN.out.gff,
        SNIPPY_RUN.out.html,
        SNIPPY_RUN.out.raw_vcf,
        SNIPPY_RUN.out.subs_vcf,
        SNIPPY_RUN.out.tab,
        SNIPPY_RUN.out.txt,
        SNIPPY_RUN.out.vcf,
        SNIPPY_RUN.out.aligned_fa_error,
        SNIPPY_RUN.out.vcf_error,
        SNIPPY_RUN.out.error
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SNIPPY_RUN.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SNIPPY_RUN.out.nf_begin,
        SNIPPY_RUN.out.nf_err,
        SNIPPY_RUN.out.nf_log,
        SNIPPY_RUN.out.nf_out,
        SNIPPY_RUN.out.nf_run,
        SNIPPY_RUN.out.nf_sh,
        SNIPPY_RUN.out.nf_trace
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SNIPPY_RUN.out.versions])
}
