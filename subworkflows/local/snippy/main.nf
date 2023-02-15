//
// snippy - Rapid variant calling from sequence reads
//
include { initOptions } from '../../../lib/nf/functions'

// reference name
reference_name = file(params.reference).getName().replaceAll(/\.fasta|\.fa|\.fna|\.gbk|\.gbff|\.gz/, '')

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
    "--maxsoft ${params.maxsoft}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
snippy_opts.subdir = "${reference_name}"

// snippy-core options
core_opts = initOptions(params.containsKey("options") ? params.options : [:], 'snippy-core')
core_opts.is_module = false
core_opts.args = "--maxhap ${params.maxhap}"
core_opts.publish_to_base = [".full.aln.gz"]
core_opts.suffix = "core-snp"
core_opts.subdir = "${reference_name}"

include { SNIPPY_RUN }  from '../../../modules/nf-core/snippy/run/main' addParams( options: snippy_opts )
include { SNIPPY_CORE }  from '../../../modules/nf-core/snippy/core/main' addParams( options: core_opts )
include { GUBBINS } from '../gubbins/main' addParams( options: [subdir:reference_name, suffix: 'core-snp', ignore: [".aln.gz"], publish_to_base: [".masked.aln.gz"]] )
include { IQTREE } from '../iqtree/main' addParams( options: [subdir:reference_name, suffix: 'core-snp', ignore: [".aln.gz"], publish_to_base: [".iqtree"]] )
include { SNPDISTS } from '../../../modules/nf-core/snpdists/main' addParams( options: [subdir:reference_name, suffix: 'core-snp.distance', publish_to_base: true] )

workflow SNIPPY {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // Run Snippy per-sample
    SNIPPY_RUN(reads, file(params.reference))
    ch_versions = ch_versions.mix(SNIPPY_RUN.out.versions.first())

    // Identify core SNPs
    SNIPPY_RUN.out.vcf.collect{meta, vcf -> vcf}.map{ vcf -> [[id:'snippy-core'], vcf]}.set{ ch_merge_vcf }
    SNIPPY_RUN.out.aligned_fa.collect{meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'snippy-core'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }
    SNIPPY_CORE(ch_snippy_core, file(params.reference))
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions.first())

    // Per-sample SNP distances
    SNPDISTS(SNIPPY_CORE.out.clean_full_aln)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    // Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(SNIPPY_CORE.out.clean_full_aln)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE.out.clean_full_aln)
        }
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }


    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
