//
// gather - Tools to gather all samples in one place
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'gather')
options.is_module = params.wf == 'gather' ? true : false
options.is_main = true
include { GATHER as GATHER_MODULE } from '../../../modules/local/bactopia/gather/main' addParams( options: options )

workflow GATHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_stats = Channel.empty()

    // Gather genomes (local, assembly, SRA/ENA)
    GATHER_MODULE(reads)
    ch_versions = ch_versions.mix(GATHER_MODULE.out.versions)

    emit:
    raw_fastq = GATHER_MODULE.out.raw_fastq
    fastq_only = GATHER_MODULE.out.fastq_only
    versions = ch_versions // channel: [ versions.yml ]
}
