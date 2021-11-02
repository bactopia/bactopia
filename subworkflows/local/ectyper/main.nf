//
// ectyper - In-silico prediction of Escherichia coli serotype
//
ectyper_opts = [
    params.verify ? "--verify" : "",
    params.print_alleles  ? "-s" : "",
    "-opid ${params.opid}",
    "-opcov ${params.opcov}",
    "-hpid ${params.hpid}",
    "-hpcov ${params.hpcov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ECTYPER as ECTYPER_MODULE } from '../../../modules/nf-core/modules/ectyper/main' addParams( options: [ args: "${ectyper_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    ECTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(ECTYPER_MODULE.out.versions.first())

    if (params.is_subworkflow) {
        ECTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ectyper'], tsv]}.set{ ch_merge_ectyper }
        CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')
    }

    emit:
    tsv = ECTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    txt = ECTYPER_MODULE.out.txt
    versions = ch_versions
}
