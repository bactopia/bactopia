//
// clonalframeml - Predict recomination events in bacterial genomes
//
clonalframeml_args = [
    "-emsim ${params.emsim}",
    "${params.clonal_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { IQTREE as START_TREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [args: "-m MFP -fast", suffix: 'start-tree', process_name: 'iqtree-cfml', is_module: "true"])
include { CLONALFRAMEML as CFML_MODULE } from '../../../modules/nf-core/modules/clonalframeml/main' addParams( options: [args: "${clonalframeml_args}", suffix: "clonalframe", is_module: "true"])

workflow CLONALFRAMEML {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_versions = Channel.empty()

    // Create a quick start tree
    START_TREE(alignment)
    ch_versions.mix(START_TREE.out.versions)

    // Run ClonalFrameML
    CFML_MODULE(START_TREE.out.aln_tree)
    ch_versions.mix(CFML_MODULE.out.versions)

    emit:
    masked_aln = CFML_MODULE.out.masked_aln
    versions = ch_versions // channel: [ versions.yml ]
}
