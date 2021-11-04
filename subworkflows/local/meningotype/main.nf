//
// meningotype - Serotyping of Neisseria meningitidis
//
meningotype_opts = [
    params.finetype ? "--finetype" : "",
    params.porB ? "--porB" : "",
    params.bast ? "--bast" : "",
    params.mlst ? "--mlst" : "",
    params.all ? "--all" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { MENINGOTYPE as MENINGOTYPE_MODULE } from '../../../modules/nf-core/modules/meningotype/main' addParams( options: [ args: "${meningotype_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow MENINGOTYPE {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    MENINGOTYPE_MODULE(fasta)
    ch_versions = ch_versions.mix(MENINGOTYPE_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        MENINGOTYPE_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'meningotype'], tsv]}.set{ ch_merge_meningotype }
        CSVTK_CONCAT(ch_merge_meningotype, 'tsv', 'tsv')
    }

    emit:
    tsv = MENINGOTYPE_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
