//
// checkm - Assess the assembly quality of your samples
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'checkm')
options.is_module = params.wf == 'checkm' ? true : false
options.args = [
    params.full_tree ? "" : "--reduced_tree",
    params.checkm_ali ? "--ali" : "",
    params.checkm_nt ? "--nt" : "",
    params.force_domain ? "--force_domain" : "",
    params.no_refinement ? "--no_refinement" : "",
    params.individual_markers ? "--individual_markers" : "",
    params.skip_adj_correction ? "--skip_adj_correction" : "",
    params.skip_pseudogene_correction ? "--skip_pseudogene_correction" : "",
    params.ignore_thresholds ? "--ignore_thresholds" : "",
    "--unique ${params.checkm_unique}",
    "--multi ${params.checkm_multi}",
    "--aai_strain ${params.aai_strain }",
    "--length ${params.checkm_length}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { CHECKM_LINEAGEWF } from '../../../modules/nf-core/modules/checkm/lineagewf/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'checkm'] )
}

workflow CHECKM {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_checkm = Channel.empty()

    CHECKM_LINEAGEWF(fasta)
    ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

    if (params.is_subworkflow) {
        CHECKM_LINEAGEWF.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'checkm'], tsv]}.set{ ch_merge_checkm }
        CSVTK_CONCAT(ch_merge_checkm, 'tsv', 'tsv')
        ch_merged_checkm = ch_merged_checkm.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    results = CHECKM_LINEAGEWF.out.results
    report = CHECKM_LINEAGEWF.out.tsv
    merged_reports = ch_merged_checkm
    versions = ch_versions // channel: [ versions.yml ]
}
