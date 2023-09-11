//
// mashtree - Quickly create a tree using Mash distances
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mashtree')
options.args = [
    "--truncLength ${params.trunclength}",
    "--sort-order ${params.sortorder}",
    "--genomesize ${params.genomesize}",
    "--mindepth ${params.mindepth}",
    "--kmerlength ${params.kmerlength}",
    "--sketch-size ${params.sketchsize}",
    params.save_sketches ? "--save-sketches sketches/" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MASHTREE as MASHTREE_MODULE } from '../../../modules/nf-core/mashtree/main' addParams( options: options )

workflow MASHTREE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    MASHTREE_MODULE(fasta)
    ch_versions = ch_versions.mix(MASHTREE_MODULE.out.versions)

    emit:
    tree = MASHTREE_MODULE.out.tree
    matrix = MASHTREE_MODULE.out.matrix
    sketches = MASHTREE_MODULE.out.sketches
    versions = ch_versions
}
