//
// phispy - Predict prophages in bacterial genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'phispy')
options.args = [
    "--number ${params.phispy_number}",
    "--min_contig_size ${params.phispy_mincontigsize}",
    "--window_size ${params.phispy_windowsize}",
    "--nonprophage_genegaps ${params.phispy_nonprophage_genegaps}",
    "--phage_genes ${params.phispy_phage_genes}",
    "--randomforest_trees ${params.phispy_randomforest_trees}",
    "${params.phispy_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PHISPY as PHISPY_MODULE } from '../../../modules/nf-core/phispy/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'phispy-concat', process_name: params.merge_folder] )

workflow PHISPY {
    take:
    gbk // channel: [ val(meta), [ gbk ] ]

    main:
    ch_versions = Channel.empty()

    PHISPY_MODULE(gbk)
    ch_versions = ch_versions.mix(PHISPY_MODULE.out.versions.first())

    // Merge results
    PHISPY_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'phispy'], tsv]}.set{ ch_merge_phispy }
    CSVTK_CONCAT(ch_merge_phispy, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = PHISPY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
