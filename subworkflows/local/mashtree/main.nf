//
// mashtree - Quickly create a tree using Mash distances
//

mashtree_args = [
    "--truncLength ${params.trunclength}",
    "--sort-order ${params.sortorder}",
    "--genomesize ${params.genomesize}",
    "--mindepth ${params.mindepth}",
    "--kmerlength ${params.kmerlength}",
    "--sketch-size ${params.sketchsize}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MASHTREE as MASHTREE_MODULE } from '../../../modules/nf-core/modules/mashtree/main' addParams( options: [ args: "${mashtree_args}", is_module: true, publish_to_base: true] )

workflow MASHTREE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    MASHTREE_MODULE(fasta)
    ch_versions = ch_versions.mix(MASHTREE_MODULE.out.versions.first())

    emit:
    tree = MASHTREE_MODULE.out.tree
    matrix = MASHTREE_MODULE.out.matrix
    versions = ch_versions
}
