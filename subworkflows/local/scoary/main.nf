//
// scoary - GWAS analysis using pangenome outputs
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'scoary')
options.is_module = params.wf == 'scoary' ? true : false
options.suffix = options.suffix ? options.suffix : 'scoary'
options.args = [
    params.permute > 0 ? "--permute ${params.permute}" : "",
    "--p_value_cutoff ${params.p_value_cutoff}",
    "--correction ${params.correction}",
    "--start_col ${params.start_col}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
SCOARY_TRAITS = params.traits ? file(params.traits) : []

include { SCOARY as SCOARY_MODULE } from '../../../modules/nf-core/modules/scoary/main' addParams( options: options )

workflow SCOARY {
    take:
    csv // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    SCOARY_MODULE(csv, SCOARY_TRAITS)
    ch_versions = ch_versions.mix(SCOARY_MODULE.out.versions)

    emit:
    csv = SCOARY_MODULE.out.csv
    versions = ch_versions
}
