//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//
if (params.use_panaroo) {
    include { PANAROO as PG_TOOL } from '../panaroo/main' addParams( options: [publish_to_base: [".aln.gz"]] )
} else if (params.use_roary) {
    include { ROARY as PG_TOOL } from '../roary/main' addParams( options: [publish_to_base: [".aln.gz"]] )
} else {
    include { PIRATE as PG_TOOL } from '../pirate/main' addParams( options: [publish_to_base: [".aln.gz"]] )
}

include { CLONALFRAMEML } from '../clonalframeml/main' addParams( options: [suffix: 'core-genome', ignore: [".aln.gz"], publish_to_base: [".masked.aln.gz"]] )
include { IQTREE as FINAL_TREE } from '../iqtree/main' addParams( options: [suffix: 'core-genome', ignore: [".aln.gz"], publish_to_base: [".iqtree"]] )
include { SNPDISTS } from '../../../modules/nf-core/snpdists/main' addParams( options: [suffix: 'core-genome.distance', publish_to_base: true] )
include { SCOARY } from '../scoary/main'

workflow PANGENOME {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ch_needs_prokka = Channel.empty()

    PG_TOOL(gff)
    ch_versions = ch_versions.mix(PG_TOOL.out.versions)

    // Per-sample SNP distances
    SNPDISTS(PG_TOOL.out.aln)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)
    
    // Identify Recombination
    if (!params.skip_recombination) {
        // Run ClonalFrameML
        CLONALFRAMEML(PG_TOOL.out.aln)
        ch_versions = ch_versions.mix(CLONALFRAMEML.out.versions)
    }

    // Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        if (params.skip_recombination) {
            FINAL_TREE(PG_TOOL.out.aln)
        } else {
            FINAL_TREE(CLONALFRAMEML.out.masked_aln)
        }
        ch_versions = ch_versions.mix(FINAL_TREE.out.versions)
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY(PG_TOOL.out.csv)
        ch_versions = ch_versions.mix(SCOARY.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
