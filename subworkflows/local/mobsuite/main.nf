//
// mobsuite - Reconstruct and annotate plasmids in bacterial assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mobsuite')
options.is_module = params.wf == 'mobsuite' ? true : false
options.args = [
    "--max_contig_size ${params.mb_max_contig_size}",
    "--min_length ${params.mb_min_contig_size}",
    "--max_plasmid_size ${params.mb_max_plasmid_size}",
    "${params.mobsuite_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MOBSUITE_RECON } from '../../../modules/nf-core/mobsuite/recon/main' addParams( options: options )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: options.is_module ? '' : 'mobsuite'] )
}

workflow MOBSUITE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_mobsuite = Channel.empty()

    MOBSUITE_RECON(fasta)
    ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions.first())

    if (params.is_subworkflow) {
        MOBSUITE_RECON.out.mobtyper_results.collect{meta, summary -> summary}.map{ summary -> [[id:'mobsuite'], summary]}.set{ ch_merge_mobsuite }
        CSVTK_CONCAT(ch_merge_mobsuite, 'tsv', 'tsv')
        ch_merged_mobsuite = ch_merged_mobsuite.mix(CSVTK_CONCAT.out.csv)
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    }

    emit:
    chromosome = MOBSUITE_RECON.out.chromosome
    contig_report = MOBSUITE_RECON.out.contig_report
    plasmids = MOBSUITE_RECON.out.plasmids
    mobtyper_results = MOBSUITE_RECON.out.mobtyper_results
    merged_reports = ch_merged_mobsuite
    versions = ch_versions
}
