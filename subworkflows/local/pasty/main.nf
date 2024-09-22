//
// pasty - in silico serogrouping of Pseudomonas aeruginosa isolates
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pasty')
options.args = [
    "--min-pident ${params.pasty_min_pident}",
    "--min-coverage ${params.pasty_min_coverage}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PASTY as PASTY_MODULE } from '../../../modules/nf-core/pasty/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'pasty-concat', process_name: params.merge_folder] )

workflow PASTY {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    PASTY_MODULE(fasta)
    ch_versions = ch_versions.mix(PASTY_MODULE.out.versions.first())

    // Merge results
    PASTY_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'pasty'], tsv]}.set{ ch_merge_pasty }
    CSVTK_CONCAT(ch_merge_pasty, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = PASTY_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PASTY_MODULE.out.blast
    versions = ch_versions
}
