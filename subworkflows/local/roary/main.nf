//
// roary - Rapid large-scale prokaryote pangenome analysis
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'roary')
options.is_module = params.wf == 'roary' ? true : false
options.args = [
    params.use_prank  ? "-e" : "-e -n",
    params.s ? "-s" : "",
    params.ap ? "-ap" : "",
    "-g ${params.g}",
    "-i ${params.i}",
    "-cd ${params.cd}",
    "-iv ${params.iv}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ROARY as ROARY_MODULE } from '../../../modules/nf-core/roary/main' addParams( options: options )

workflow ROARY {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ROARY_MODULE(gff)
    ch_versions = ch_versions.mix(ROARY_MODULE.out.versions)

    emit:
    aln = ROARY_MODULE.out.aln
    csv = ROARY_MODULE.out.csv
    versions = ch_versions
}
