//
// pirate - Pangenome toolbox for bacterial genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pirate')
options.is_module = params.wf == 'pirate' ? true : false
options.args = [
    params.para_off ? "--para_off" : "",
    params.pan_opt ? "--pan-opt '${params.pan_opt}'" : "",
    params.z ? "-z 2" : "",
    "--steps ${params.steps}",
    "--features ${params.features}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PIRATE as PIRATE_MODULE } from '../../../modules/nf-core/modules/pirate/main' addParams( options: options )

workflow PIRATE {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    PIRATE_MODULE(gff)
    ch_versions = ch_versions.mix(PIRATE_MODULE.out.versions.first())

    emit:
    aln = PIRATE_MODULE.out.aln
    csv = PIRATE_MODULE.out.csv
    versions = ch_versions
}
