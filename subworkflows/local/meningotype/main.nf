//
// meningotype - Serotyping of Neisseria meningitidis
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'meningotype')
options.args = [
    params.finetype ? "--finetype" : "",
    params.porB ? "--porB" : "",
    params.bast ? "--bast" : "",
    params.mlst ? "--mlst" : "",
    params.all ? "--all" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../../modules/nf-core/meningotype/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'meningotype'] )

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MENINGOTYPE_MODULE(fasta)
    ch_versions = ch_versions.mix(MENINGOTYPE_MODULE.out.versions.first())

    // Merge results
    MENINGOTYPE_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'meningotype'], tsv]}.set{ ch_merge_meningotype }
    CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
