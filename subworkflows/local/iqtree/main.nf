//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'iqtree')
options.is_module = params.wf == 'iqtree' ? true : false
options.args = [
    params.asr ? "-asr" : "",
    "-m ${params.iqtree_model}",
    "-bb ${params.bb}",
    "-alrt ${params.alrt}",
    "-wbt -wbtl -alninfo ${params.iqtree_opts}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { IQTREE as IQTREE_MODULE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: options )

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
