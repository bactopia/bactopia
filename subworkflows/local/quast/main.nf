//
// quast - A module for assessing the quality of assembled contigs
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'quast')
options.args = [
    "--contig-thresholds '${params.contig_thresholds}'",
    "--plots-format ${params.plots_format}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { QUAST as QUAST_MODULE } from '../../../modules/nf-core/quast/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'seroba-concat', process_name: params.merge_folder] )

workflow QUAST {
    take:
    fasta // channel: [ val(meta), [ fasta ], [ meta_files ] ]

    main:
    ch_versions = Channel.empty()

    QUAST_MODULE(fasta)
    ch_versions = ch_versions.mix(QUAST_MODULE.out.versions.first())

    // Merge results
    QUAST_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'quast'], tsv]}.set{ ch_merge_quast}
    CSVTK_CONCAT(ch_merge_quast, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = QUAST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
