//
// busco - Assembly completeness based on evolutionarily informed expectations
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'busco')
options.is_module = params.wf == 'busco' ? true : false
options.args = [
    params.metaeuk_parameters ? "--metaeuk_parameters '${params.metaeuk_parameters}'" : "",
    params.metaeuk_rerun_parameters ? "--metaeuk_rerun_parameters '${params.metaeuk_rerun_parameters}'" : "",
    params.use_augustus ? "--augustus" : "",
    params.augustus_parameters ? "--augustus_parameters '${params.augustus_parameters}'" : "",
    params.augustus_species ? "--augustus_species '${params.augustus_species}'" : "",
    params.augustus_long ? "--long" : "",
    "--evalue ${params.busco_evalue}",
    "--limit ${params.busco_limit}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { BUSCO as BUSCO_MODULE } from '../../../modules/nf-core/busco/main' addParams( options: options )

if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'busco'] )
}

workflow BUSCO {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_lineages = Channel.from(params.busco_lineage.split(','))
    ch_merged_busco = Channel.empty()

    BUSCO_MODULE(fasta, ch_lineages)
    ch_versions = ch_versions.mix(BUSCO_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        BUSCO_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'busco'], tsv]}.set{ ch_merge_busco }
        CSVTK_CONCAT(ch_merge_busco, 'tsv', 'tsv')
        ch_merged_busco = ch_merged_busco.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    results = BUSCO_MODULE.out.results
    tsv = BUSCO_MODULE.out.tsv
    merged_tsv = ch_merged_busco
    versions = ch_versions // channel: [ versions.yml ]
}
