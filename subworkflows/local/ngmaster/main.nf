//
// ngmaster - Multi-antigen sequence typing for Neisseria gonorrhoeae
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ngmaster')
options.args = [
    params.csv ? "--csv" : ""
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { NGMASTER as NGMASTER_MODULE } from '../../../modules/nf-core/ngmaster/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'ngmaster'] )

workflow NGMASTER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    NGMASTER_MODULE(fasta)
    ch_versions = ch_versions.mix(NGMASTER_MODULE.out.versions.first())

    // Merge results
    NGMASTER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ngmaster'], tsv]}.set{ ch_merge_ngmaster }
    CSVTK_CONCAT(ch_merge_ngmaster, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = NGMASTER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions
}
