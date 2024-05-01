//
// mlst - Automatic MLST calling from assembled contigs
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'mlst')
options.args = [
    params.nopath ? "--nopath" : "",
    params.scheme ? "--scheme ${params.scheme}" : "",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}",
    "--minscore ${params.minscore}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MLST as MLST_MODULE } from '../../../modules/nf-core/mlst/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [args: '--no-header-row', logs_subdir: 'mlst-concat', process_name: params.merge_folder] )

workflow MLST {
    take:
    fasta // channel: [ val(meta), [ reads ] ]
    db // channel: [ mlst_db ]

    main:
    ch_versions = Channel.empty()

    if (params.mlst_db) {
        MLST_MODULE(fasta, file(params.mlst_db))
    } else {
        MLST_MODULE(fasta, db)
    }

    ch_versions = ch_versions.mix(MLST_MODULE.out.versions.first())

    // Merge results
    MLST_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'mlst'], tsv]}.set{ ch_merge_mlst }
    CSVTK_CONCAT(ch_merge_mlst, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = MLST_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
