//
// ectyper - In-silico prediction of Escherichia coli serotype
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ectyper')
options.args = [
    params.verify ? "--verify" : "",
    params.print_alleles  ? "-s" : "",
    "-opid ${params.opid}",
    "-opcov ${params.opcov}",
    "-hpid ${params.hpid}",
    "-hpcov ${params.hpcov}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ECTYPER as ECTYPER_MODULE } from '../../../modules/nf-core/ectyper/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'ectyper-concat', process_name: params.merge_folder] )

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    ECTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(ECTYPER_MODULE.out.versions.first())

    // Merge results
    ECTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ectyper'], tsv]}.set{ ch_merge_ectyper }
    CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = ECTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    txt = ECTYPER_MODULE.out.txt
    versions = ch_versions
}
