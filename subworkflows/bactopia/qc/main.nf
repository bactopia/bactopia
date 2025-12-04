//
// qc - Quality control of Illumina and ONT reads
//
nextflow.preview.types = true

include { QC as QC_MODULE } from '../../../modules/bactopia/qc/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow QC {
    take:
    reads: Channel<Tuple<Map, Set<Path>, Set<Path>>>
    adapters: Path?
    phix: Path?

    main:
    QC_MODULE(reads, adapters, phix)

    emit:
    // Individual outputs
    fastq: Channel<Tuple<Map, Set<Path>, Set<Path>>> = QC_MODULE.out.fastq
    fastq_only: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.fastq_only
    txt: Channel<Tuple<Map, Path>> = QC_MODULE.out.txt
    error: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error
    error_fastq: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error_fastq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        QC_MODULE.out.fastq_only,
        QC_MODULE.out.txt,
        QC_MODULE.out.supplemental,
        QC_MODULE.out.error,
        QC_MODULE.out.error_fastq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.versions])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.logs])
}
