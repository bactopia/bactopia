//
// panaroo - Pipeline for pangenome investigations
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'panaroo')
options.is_module = params.wf == 'panaroo' ? true : false
options.args = [
    params.merge_paralogs ? "--merge_paralogs" : "",
    "--clean-mode ${params.panaroo_mode}",
    "--threshold ${params.panaroo_threshold}",
    "--family_threshold ${params.panaroo_family_threshold}",
    "--len_dif_percent ${params.len_dif_percent}",
    "--alignment ${params.panaroo_alignment}",
    "--aligner ${params.panaroo_aligner}",
    "--core_threshold ${params.panaroo_core_threshold}",
    params.panaroo_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PANAROO_RUN } from '../../../modules/nf-core/panaroo/run/main' addParams( options: options )

workflow PANAROO {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    PANAROO_RUN(gff)
    ch_versions = ch_versions.mix(PANAROO_RUN.out.versions)

    emit:
    aln = PANAROO_RUN.out.aln
    csv = PANAROO_RUN.out.csv
    panaroo_csv = PANAROO_RUN.out.panaroo_csv
    versions = ch_versions
}
