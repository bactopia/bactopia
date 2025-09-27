//
// snippy - Rapid variant calling from sequence reads
//
include { SNIPPY_RUN  } from '../../modules/snippy/run/main'
include { SNIPPY_CORE } from '../../modules/snippy/core/main'
include { GUBBINS     } from '../gubbins/main'
include { IQTREE      } from '../iqtree/main'
include { SNPDISTS as SNPDISTS_UNMASKED } from '../../modules/snpdists/main'
include { SNPDISTS as SNPDISTS_MASKED   } from '../../modules/snpdists/main'
workflow SNIPPY {
    take:
    reads     // channel: [ val(meta), [ reads ] ]
    reference // channel: [ val(meta), [ fasta ] ]
    mask

    main:
    ch_versions = Channel.empty()

    // Run Snippy per-sample
    SNIPPY_RUN(reads, reference)
    ch_versions = ch_versions.mix(SNIPPY_RUN.out.versions)

    emit:
    aligned_fa               = SNIPPY_RUN.out.aligned_fa
    annotated_vcf            = SNIPPY_RUN.out.annotated_vcf
    bam                      = SNIPPY_RUN.out.bam
    bai                      = SNIPPY_RUN.out.bai
    bed                      = SNIPPY_RUN.out.bed
    consensus_fa             = SNIPPY_RUN.out.consensus_fa
    consensus_subs_fa        = SNIPPY_RUN.out.consensus_subs_fa
    consensus_subs_masked_fa = SNIPPY_RUN.out.consensus_subs_masked_fa
    coverage                 = SNIPPY_RUN.out.coverage
    csv                      = SNIPPY_RUN.out.csv
    filt_vcf                 = SNIPPY_RUN.out.filt_vcf
    gff                      = SNIPPY_RUN.out.gff
    html                     = SNIPPY_RUN.out.html
    raw_vcf                  = SNIPPY_RUN.out.raw_vcf
    subs_vcf                 = SNIPPY_RUN.out.subs_vcf
    tab                      = SNIPPY_RUN.out.tab
    txt                      = SNIPPY_RUN.out.txt
    vcf                      = SNIPPY_RUN.out.vcf
    logs                     = SNIPPY_RUN.out.logs
    nf_logs                  = SNIPPY_RUN.out.nf_begin.mix(SNIPPY_RUN.out.nf_err, SNIPPY_RUN.out.nf_log, SNIPPY_RUN.out.nf_out, SNIPPY_RUN.out.nf_run, SNIPPY_RUN.out.nf_sh, SNIPPY_RUN.out.nf_trace)
    versions                 = ch_versions
}
