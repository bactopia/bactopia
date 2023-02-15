//
// qc - Quality control of Illumina and ONT reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'qc')
options.is_module = params.wf == 'qc' ? true : false
options.is_main = true
ADAPTERS = params.adapters ? file(params.adapters) : []
PHIX = params.phix ? file(params.phix) : []

include { QC as QC_MODULE } from '../../../modules/local/bactopia/qc/main' addParams( options: options )

workflow QC {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_stats = Channel.empty()

    // QC reads
    QC_MODULE(reads, ADAPTERS, PHIX)
    ch_versions = ch_versions.mix(QC_MODULE.out.versions)

    emit:
    fastq = QC_MODULE.out.fastq
    fastq_assembly = QC_MODULE.out.fastq_assembly
    versions = ch_versions // channel: [ versions.yml ]
}
