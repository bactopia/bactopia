//
// qc - Quality control of Illumina and ONT reads
//
include { QC as QC_MODULE } from '../../modules/bactopia/qc/main'

workflow QC {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Set up input files
    ADAPTERS = params.adapters ? file(params.adapters) : []
    PHIX = params.phix ? file(params.phix) : []

    // QC reads
    QC_MODULE(reads, ADAPTERS, PHIX)
    ch_versions = ch_versions.mix(QC_MODULE.out.versions)
    ch_logs = ch_logs.mix(QC_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QC_MODULE.out.nf_logs)

    emit:
    fastq = QC_MODULE.out.fastq
    fastq_only = QC_MODULE.out.fastq_only
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
