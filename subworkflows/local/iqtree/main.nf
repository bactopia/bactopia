//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
iqtree_args = [
    params.asr ? "-asr" : "",
    "-m ${params.m}",
    "-bb ${params.bb}",
    "-alrt ${params.alrt}",
    "-wbt -wbtl -alninfo ${params.iqtree_opts}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { IQTREE as IQTREE_MODULE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [ args: "${iqtree_args}", is_module: true] )

workflow IQTREE {
    take:
    aln // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    IQTREE_MODULE(aln)
    ch_versions = ch_versions.mix(IQTREE_MODULE.out.versions.first())

    emit:
    phylogeny = IQTREE_MODULE.out.phylogeny
    versions = ch_versions // channel: [ versions.yml ]
}
