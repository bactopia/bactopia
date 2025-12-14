/**
 * Call variants against a reference genome using Snippy.
 *
 * This subworkflow performs rapid haploid variant calling from bacterial sequence reads
 * using [Snippy](https://github.com/tseemann/snippy). It maps reads to a reference genome,
 * identifies SNPs and indels, and generates consensus sequences. The tool produces multiple
 * output formats including VCF, aligned FASTA, and annotated variants for downstream
 * phylogenetic analysis with snippy-core.
 *
 * @status stable
 * @keywords variant calling, snp, reference mapping, phylogenetics, outbreak
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation snippy
 *
 * @modules snippy_run
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Sequencing reads (paired-end or single-end FASTQ) for variant calling
 *
 * @input reference
 * Reference genome in GenBank format (preferred, for annotation) or FASTA format
 *
 * @output aligned_fa               Per-sample aligned FASTA for core-genome analysis
 * @output vcf                      Filtered variant calls in VCF format
 * @output aligned_fa_error         Aligned FASTA from samples that encountered errors
 * @output vcf_error                VCF from samples that encountered errors
 * @output error                    Error messages from failed samples
 * @output annotated_vcf            Annotated variants with gene/product information
 * @output bam                      Aligned reads in BAM format
 * @output bai                      BAM index file
 * @output bed                      Genomic positions of variants in BED format
 * @output consensus_fa             Consensus sequence in FASTA format
 * @output consensus_subs_fa        Consensus with only substitutions applied
 * @output consensus_subs_masked_fa Consensus with low-coverage regions masked as N
 * @output coverage                 Per-base coverage statistics
 * @output csv                      Variants in CSV format for easy parsing
 * @output filt_vcf                 Quality-filtered VCF file
 * @output gff                      Annotated variants in GFF3 format
 * @output html                     Summary report in HTML format
 * @output raw_vcf                  Unfiltered variant calls
 * @output subs_vcf                 Substitution-only variants in VCF format
 * @output tab                      Tab-delimited variant summary
 * @output txt                      Text summary of variant calling statistics
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
    aligned_fa: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.aligned_fa
    vcf: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.vcf
    aligned_fa_error: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.aligned_fa_error
    vcf_error: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.vcf_error
    error: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.error
    annotated_vcf: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.annotated_vcf
    bam: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.bam
    bai: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.bai
    bed: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.bed
    consensus_fa: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.consensus_fa
    consensus_subs_fa: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.consensus_subs_fa
    consensus_subs_masked_fa: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.consensus_subs_masked_fa
    coverage: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.coverage
    csv: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.csv
    filt_vcf: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.filt_vcf
    gff: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.gff
    html: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.html
    raw_vcf: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.raw_vcf
    subs_vcf: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.subs_vcf
    tab: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.tab
    txt: Channel<Tuple<Map, Set<Path>>> = SNIPPY_RUN.out.txt

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
