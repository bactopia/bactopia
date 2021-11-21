//
// roary - Rapid large-scale prokaryote pangenome analysis
//
roary_args = [
    params.use_prank  ? "-e" : "-e -n",
    params.s ? "-s" : "",
    params.ap ? "-ap" : "",
    "-g ${params.g}",
    "-i ${params.i}",
    "-cd ${params.cd}",
    "-iv ${params.iv}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ROARY as ROARY_MODULE } from '../../../modules/nf-core/modules/roary/main' addParams( options: [ args: "${roary_args}", is_module: true] )

workflow ROARY {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ROARY_MODULE(gff)
    ch_versions = ch_versions.mix(ROARY_MODULE.out.versions.first())

    emit:
    aln = ROARY_MODULE.out.aln
    csv = ROARY_MODULE.out.csv
    versions = ch_versions
}
