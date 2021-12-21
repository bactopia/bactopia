//
// clonalframeml - Predict recomination events in bacterial genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'clonalframeml')
options.args = ["-emsim ${params.emsim}", "${params.clonal_opts}"].join(' ').replaceAll("\\s{2,}", " ").trim()
options.is_module = params.wf == 'clonalframeml' ? true : false
options.suffix = options.suffix ? options.suffix : 'clonalframe'

include { IQTREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [args: "-m MFP -fast", suffix: 'start-tree', process_name: 'iqtree-fast', is_module: options.is_module, ignore: options.ignore])
include { CLONALFRAMEML as CLONALFRAME } from '../../../modules/nf-core/modules/clonalframeml/main' addParams( options: options)

workflow CLONALFRAMEML {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_versions = Channel.empty()

    // Create a quick start tree
    IQTREE(alignment)
    ch_versions = ch_versions.mix(IQTREE.out.versions)

    // Run ClonalFrameML
    CLONALFRAME(IQTREE.out.aln_tree)
    ch_versions = ch_versions.mix(CLONALFRAME.out.versions)

    emit:
    masked_aln = CLONALFRAME.out.masked_aln
    versions = ch_versions // channel: [ versions.yml ]
}
