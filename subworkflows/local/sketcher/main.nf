//
// sketcher - Assortment of tools for sketching sequences
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'sketcher')
options.is_module = params.wf == 'sketcher' ? true : false
options.is_main = true

// args -> Mash
options.args = [
    "-s ${params.sketch_size}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

// args2 -> Sourmash
options.args2 = [
    "-p k=21,k=31,k=51,abund,scaled=${params.sourmash_scale}"
].join(' ').replaceAll("\\s{2,}", " ").trim()


include { SKETCHER as SKETCHER_MODULE } from '../../../modules/local/bactopia/sketcher/main' addParams( options: options )

workflow SKETCHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_stats = Channel.empty()

    // Sketch FASTQs
    SKETCHER_MODULE(reads)
    ch_versions = ch_versions.mix(SKETCHER_MODULE.out.versions)

    emit:
    sig = SKETCHER_MODULE.out.sig
    versions = ch_versions // channel: [ versions.yml ]
}
