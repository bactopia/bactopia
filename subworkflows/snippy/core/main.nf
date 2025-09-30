//
// snippy - Rapid variant calling from sequence reads
//
include { SNIPPY_CORE as SNIPPY_CORE_MODULE  } from '../../../modules/snippy/core/main'
include { SNPDISTS                           } from '../../snpdists/main'
workflow SNIPPY_CORE {
    take:
    alignments // channel: [ val(meta), [ reads ] ]
    reference  // channel: [ val(meta), [ fasta ] ]
    mask

    main:
    SNIPPY_CORE_MODULE(alignments, reference, mask)

    // Per-sample SNP distances
    SNIPPY_CORE_MODULE.out.clean_full_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-snp.distance", process_name: "snpdists"], aln]}.set{ ch_unmasked_aln }
    SNPDISTS(ch_unmasked_aln)

    emit:
    aln = SNIPPY_CORE_MODULE.out.aln
    full_aln = SNIPPY_CORE_MODULE.out.full_aln
    clean_full_aln = SNIPPY_CORE_MODULE.out.clean_full_aln
    tab = SNIPPY_CORE_MODULE.out.tab
    vcf = SNIPPY_CORE_MODULE.out.vcf
    txt = SNIPPY_CORE_MODULE.out.txt
    samples = SNIPPY_CORE_MODULE.out.samples
    tsv = SNPDISTS.out.tsv

    // Generic aggregate outputs
    results = SNIPPY_CORE_MODULE.out.supplemental.mix(
        SNIPPY_CORE_MODULE.out.aln,
        SNIPPY_CORE_MODULE.out.full_aln,
        SNIPPY_CORE_MODULE.out.clean_full_aln,
        SNIPPY_CORE_MODULE.out.tab,
        SNIPPY_CORE_MODULE.out.vcf,
        SNIPPY_CORE_MODULE.out.txt,
        SNIPPY_CORE_MODULE.out.samples,
        SNPDISTS.out.results
    )
    logs = SNIPPY_CORE_MODULE.out.logs.mix(
        SNPDISTS.out.logs
    )
    nf_logs = SNIPPY_CORE_MODULE.out.nf_begin.mix(
        SNIPPY_CORE_MODULE.out.nf_err,
        SNIPPY_CORE_MODULE.out.nf_log,
        SNIPPY_CORE_MODULE.out.nf_out,
        SNIPPY_CORE_MODULE.out.nf_run,
        SNIPPY_CORE_MODULE.out.nf_sh,
        SNIPPY_CORE_MODULE.out.nf_trace,
        SNPDISTS.out.nf_logs
    )
    versions = SNIPPY_CORE_MODULE.out.versions.mix(
        SNPDISTS.out.versions
    )
}
