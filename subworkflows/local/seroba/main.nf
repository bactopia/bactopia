//
// seroba - Serotyping of Streptococcus pneumoniae from sequence reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'seroba')
options.args = [
    "--coverage ${params.seroba_coverage}",
    params.seroba_noclean ? "--noclean" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { SEROBA_RUN } from '../../../modules/nf-core/seroba/run/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'seroba-concat', process_name: params.merge_folder] )

workflow SEROBA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    SEROBA_RUN(fasta)
    ch_versions = ch_versions.mix(SEROBA_RUN.out.versions.first())

    // Merge results
    SEROBA_RUN.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'seroba'], tsv]}.set{ ch_merge_seroba }
    CSVTK_CONCAT(ch_merge_seroba, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SEROBA_RUN.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
