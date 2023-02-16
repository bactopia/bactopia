//
// seqsero2 - Salmonella serotype prediction from reads or assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'seqsero2')
options.args = [
    "-m ${params.run_mode}",
    params.input_type == "assembly" ? "-t 4" : "-t 2",
    "-b ${params.bwa_mode}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
include { SEQSERO2 as SEQSERO2_MODULE } from '../../../modules/nf-core/seqsero2/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'seqsero2'] )

workflow SEQSERO2 {
    take:
    seqs // channel: [ val(meta), [ fastqs or assemblies ] ]

    main:
    ch_versions = Channel.empty()

    SEQSERO2_MODULE(seqs)
    ch_versions = ch_versions.mix(SEQSERO2_MODULE.out.versions.first())

    // Merge results
    SEQSERO2_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'seqsero2'], tsv]}.set{ ch_merge_seqsero2 }
    CSVTK_CONCAT(ch_merge_seqsero2, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = SEQSERO2_MODULE.out.tsv
    txt = SEQSERO2_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
