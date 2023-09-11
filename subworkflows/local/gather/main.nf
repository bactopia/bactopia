//
// gather - Tools to gather all samples in one place
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'gather')
include { GATHER as GATHER_MODULE } from '../../../modules/local/bactopia/gather/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'meta-concat', process_name: params.merge_folder] )

workflow GATHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // Gather genomes (local, assembly, SRA/ENA)
    GATHER_MODULE(reads)
    ch_versions = ch_versions.mix(GATHER_MODULE.out.versions)

    // Merge meta values for each sample
    GATHER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'meta'], tsv]}.set{ ch_merge_stats }
    CSVTK_CONCAT(ch_merge_stats, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    raw_fastq = GATHER_MODULE.out.raw_fastq
    fastq_only = GATHER_MODULE.out.fastq_only
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
