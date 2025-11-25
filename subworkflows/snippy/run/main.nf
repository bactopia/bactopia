//
// snippy - Rapid variant calling from sequence reads
//
nextflow.preview.types = true

include { SNIPPY_RUN  } from '../../../modules/snippy/run/main'
workflow SNIPPY {
    take:
    reads     : Channel<Tuple<Map, List<Path>>>
    reference : Tuple<Map, Path>

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
    results: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.aligned_fa.mix(
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
    )
    logs: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.nf_begin.mix(
        SNIPPY_RUN.out.nf_err,
        SNIPPY_RUN.out.nf_log,
        SNIPPY_RUN.out.nf_out,
        SNIPPY_RUN.out.nf_run,
        SNIPPY_RUN.out.nf_sh,
        SNIPPY_RUN.out.nf_trace
    )
    versions: Channel<Tuple<Map, Path>> = SNIPPY_RUN.out.versions
}
