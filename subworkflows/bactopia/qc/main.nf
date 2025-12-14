/**
 * Perform comprehensive quality control on sequencing reads.
 *
 * This subworkflow processes raw sequencing reads through a comprehensive quality control pipeline.
 * It adapts to different read types (Illumina paired-end, single-end, and Nanopore long reads),
 * performing adapter/PhiX removal, error correction, quality filtering, and coverage reduction.
 * Generates detailed quality reports using [FastQC](https://github.com/s-andrews/FastQC) and [NanoPlot](https://github.com/wdecoster/NanoPlot).
 *
 * Uses explicit positional tuple slots for reads:
 * - Input/Output: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords quality control, adapters, error correction, subsampling, fastq, illumina, nanopore
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, conditional-logic, path-workarounds
 * @citation bbmap, fastp, fastqc, fastq-scan, lighter, nanoplot, nanoq, porechop, rasusa
 *
 * @modules qc
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) or original assembly (if simulated)
 *
 * @input adapters
 * Optional adapter sequences in FASTA format for removal from Illumina reads
 *
 * @input phix
 * Optional PhiX sequences in FASTA format for removal from Illumina reads
 *
 * @output reads        Tuple with explicit read slots: (meta, r1, r2, se, lr) where each is Path?
 * @output error_fastq  Reads preserved from samples that failed QC for debugging
 * @output supplemental QC reports, quality metrics, and original/final FASTQ comparisons
 * @output error        Error messages from QC failures
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { QC as QC_MODULE } from '../../../modules/bactopia/qc/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow QC {
    take:
    samples: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    adapters: Path?
    phix: Path?

    main:
    QC_MODULE(samples, adapters, phix)

    emit:
    // Individual outputs
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>> = QC_MODULE.out.reads
    error: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error
    error_fastq: Channel<Tuple<Map, Set<Path>>> = QC_MODULE.out.error_fastq

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        QC_MODULE.out.supplemental,
        QC_MODULE.out.error,
        QC_MODULE.out.error_fastq
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.versions])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([QC_MODULE.out.logs])
}
