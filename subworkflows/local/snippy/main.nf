//
// snippy - Rapid variant calling from sequence reads
//
include { initOptions } from '../../../lib/nf/functions'

// snippy options
snippy_opts = initOptions(params.containsKey("options") ? params.options : [:], 'snippy')
snippy_opts.is_module = params.wf == 'snippy' ? true : false
snippy_opts.args = [
    params.bwaopt ? "--bwaopt ${params.bwaopt}" : "",
    params.fbopt ? "--fbopt ${params.fbopt}" : "",
    "--mapqual ${params.mapqual}",
    "--basequal ${params.basequal}",
    "--mincov ${params.mincov}",
    "--minfrac ${params.minfrac}",
    "--minqual ${params.minqual}",
    "--maxsoft ${params.maxsoft}",
    params.snippy_opts ? "${params.snippy_opts}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
snippy_opts.subdir = params.run_name
snippy_opts.logs_use_prefix = true

// snippy-core options
MASK = params.mask ? file(params.mask) : []
core_opts = initOptions(params.containsKey("options") ? params.options : [:], 'snippy-core')
core_opts.is_module = false
core_opts.args = [
    "--maxhap ${params.maxhap}",
    "--mask-char ${params.mask_char}",
    params.snippy_core_opts ? "${params.snippy_core_opts}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
core_opts.publish_to_base = [".full.aln.gz", ".samples.txt"]
core_opts.suffix = "core-snp"

include { SNIPPY_RUN as SNIPPY_RUN_MODULE }  from '../../../modules/nf-core/snippy/run/main' addParams( options: snippy_opts )
include { SNIPPY_CORE as SNIPPY_CORE_MODULE }  from '../../../modules/nf-core/snippy/core/main' addParams( options: core_opts )
include { GUBBINS } from '../gubbins/main' addParams( options: [suffix: 'core-snp', ignore: [".aln.gz"]] )
include { IQTREE } from '../iqtree/main' addParams( options: [suffix: 'core-snp', ignore: [".aln.gz"]] )
include { SNPDISTS as SNPDISTS_UNMASKED } from '../../../modules/nf-core/snpdists/main' addParams( options: [suffix: 'core-snp.distance', publish_to_base: true, logs_subdir: "unmasked-aln"] )
include { SNPDISTS as SNPDISTS_MASKED } from '../../../modules/nf-core/snpdists/main' addParams( options: [suffix: 'core-snp.masked.distance', publish_to_base: true, logs_subdir: "masked-aln"] )

workflow SNIPPY {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    reference // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run Snippy per-sample
    SNIPPY_RUN_MODULE(reads, reference)
    ch_versions = ch_versions.mix(SNIPPY_RUN_MODULE.out.versions.first())

    // Identify core SNPs
    SNIPPY_RUN_MODULE.out.vcf.collect{meta, vcf -> vcf}.map{ vcf -> [[id:'snippy-core'], vcf]}.set{ ch_merge_vcf }
    SNIPPY_RUN_MODULE.out.aligned_fa.collect{meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'snippy-core'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }
    SNIPPY_CORE_MODULE(ch_snippy_core, reference, MASK)
    ch_versions = ch_versions.mix(SNIPPY_CORE_MODULE.out.versions.first())

    // Per-sample SNP distances
    SNPDISTS_UNMASKED(SNIPPY_CORE_MODULE.out.clean_full_aln)
    ch_versions = ch_versions.mix(SNPDISTS_UNMASKED.out.versions)

    // Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(SNIPPY_CORE_MODULE.out.clean_full_aln)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
        SNPDISTS_MASKED(GUBBINS.out.masked_aln)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE_MODULE.out.clean_full_aln)
        }
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

workflow SNIPPY_RUN {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // Run Snippy per-sample
    SNIPPY_RUN_MODULE(reads, file(params.reference))
    ch_versions = ch_versions.mix(SNIPPY_RUN_MODULE.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

workflow SNIPPY_CORE {
    take:
    snippy_run // channel: [ val(meta), [ vcf ], [ aligned_fa ] ]

    main:
    ch_versions = Channel.empty()

    // Identify core SNPs
    snippy_run.vcf.collect{meta, vcf -> vcf}.map{ vcf -> [[id:'snippy-core'], vcf]}.set{ ch_merge_vcf }
    snippy_run.aligned_fa.collect{meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'snippy-core'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }
    SNIPPY_CORE_MODULE(ch_snippy_core, file(params.reference), MASK)
    ch_versions = ch_versions.mix(SNIPPY_CORE_MODULE.out.versions.first())

    // Per-sample SNP distances
    SNPDISTS_UNMASKED(SNIPPY_CORE_MODULE.out.clean_full_aln)
    ch_versions = ch_versions.mix(SNPDISTS_UNMASKED.out.versions)

    // Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(SNIPPY_CORE_MODULE.out.clean_full_aln)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
        SNPDISTS_MASKED(GUBBINS.out.masked_aln)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE_MODULE.out.clean_full_aln)
        }
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
