//
// ngmaster - Multi-antigen sequence typing for Neisseria gonorrhoeae
//
ngmaster_opts = [
    params.csv ? "--csv" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { NGMASTER as NGMASTER_MODULE } from '../../../modules/nf-core/modules/ngmaster/main' addParams( options: [ args: "${ngmaster_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow NGMASTER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    NGMASTER_MODULE(fasta)
    ch_versions = ch_versions.mix(NGMASTER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        NGMASTER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ngmaster'], tsv]}.set{ ch_merge_ngmaster }
        CSVTK_CONCAT(ch_merge_ngmaster, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }

    emit:
    tsv = NGMASTER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
