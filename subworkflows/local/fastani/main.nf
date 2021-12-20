//
// fastani - fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
//
fastani_opts = [
    "--kmer ${params.kmer}",
    "--fragLen ${params.frag_len}",
    "--minFraction ${params.min_fraction}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { FASTANI as FASTANI_MODULE } from '../../../modules/nf-core/modules/fastani/main' addParams( options: [ args: "${fastani_opts}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow FASTANI {
    take:
    query // channel: [ val(meta), [ fasta ] ]
    reference // channel: [ fastas ]

    main:
    ch_versions = Channel.empty()
    ch_fastani_reference = Channel.empty()
    query.collect{meta, fasta -> fasta}.map{ fasta -> [[id:'query'], fasta]}.set{ ch_fastani_query }
    if (params.skip_pairwise) {
        // All against each reference
        ch_fastani_reference = reference.map{meta, fasta -> fasta}
    } else {
        // All by All
        ch_fastani_reference = query.map{meta, fasta -> fasta}
    }

    FASTANI_MODULE ( ch_fastani_query, ch_fastani_reference )
    if (params.is_subworkflow) {
        FASTANI_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'fastani'], tsv]}.set{ ch_merge_fastani }
        CSVTK_CONCAT(ch_merge_fastani, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }

    ch_versions = ch_versions.mix(FASTANI_MODULE.out.versions.first())

    emit:
    tsv = FASTANI_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
