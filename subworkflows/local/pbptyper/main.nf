//
// pbptyper - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pbptyper')
options.args = [
    "--min_pident ${params.pbptyper_min_pident}",
    "--min_coverage ${params.pbptyper_min_coverage}",
    "--min_ani ${params.pbptyper_min_ani}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PBPTYPER as PBPTYPER_MODULE } from '../../../modules/nf-core/pbptyper/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [process_name: 'pbptyper'] )

workflow PBPTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    PBPTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(PBPTYPER_MODULE.out.versions.first())

    // Merge results
    PBPTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'pbptyper'], tsv]}.set{ ch_merge_pbptyper }
    CSVTK_CONCAT(ch_merge_pbptyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = PBPTYPER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    blast = PBPTYPER_MODULE.out.blast
    versions = ch_versions
}
