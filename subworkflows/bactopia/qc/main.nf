/**
 * Perform comprehensive quality control on sequencing reads.
 *
 * This subworkflow processes raw sequencing reads through a comprehensive quality control pipeline.
 * It adapts to different read types:
 * - **Illumina:** Adapter/PhiX removal ([Fastp](https://github.com/OpenGene/fastp) or
 *   [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)), Error Correction
 *   ([Lighter](https://github.com/mourisl/Lighter)), and Subsampling ([Rasusa](https://github.com/mbhall88/rasusa))
 * - **Nanopore:** Adapter removal ([Porechop](https://github.com/rrwick/Porechop)), Quality filtering
 *   ([Nanoq](https://github.com/esteinig/nanoq)), and Subsampling ([Rasusa](https://github.com/mbhall88/rasusa))
 * - **Hybrid:** Processes both short and long reads through their respective pipelines
 * - **Assembly:** Passes through simulated reads from assemblies
 *
 * Generates quality metrics using [fastq-scan](https://github.com/rpetit3/fastq-scan) and optional
 * quality reports using [FastQC](https://github.com/s-andrews/FastQC) (Illumina) and
 * [NanoPlot](https://github.com/wdecoster/NanoPlot) (ONT).
 *
 * @status stable
 * @keywords quality control, adapters, error correction, subsampling, fastq, illumina, nanopore, fastp, bbduk, nanoq
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation,conditional-input,conditional-logic
 * @citation bbtools, fastp, fastqc, fastq_scan, lighter, nanoplot, nanoq, porechop, rasusa
 *
 * @modules bactopia_qc
 *
 * @input record(meta, r1?, r2?, se?, lr?, assembly?)
 * - `meta`: Groovy Record containing sample information (must include `runtype`, `genome_size`, `species`)
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT)
 * - `assembly?`: Assembly file (FASTA) for assembly-based simulations
 *
 * @input adapters
 * Optional adapter sequences in FASTA format for removal from Illumina reads
 *
 * @input phix
 * Optional PhiX sequences in FASTA format for removal from Illumina reads
 *
 * @output sample_outputs
 * - `reads_grouped`: All output FASTQs for publishing
 * - `supplemental`: QC reports (FastQC/NanoPlot), JSON metrics, and error FASTQs if QC failed
 * - `error`: Captured error messages if QC failed (e.g., reads empty after trimming)
 *
 * @output run_outputs
 *
 * @output reads
 * - `r1?`: QC-filtered Illumina R1 reads
 * - `r2?`: QC-filtered Illumina R2 reads
 * - `se?`: QC-filtered single-end reads
 * - `lr?`: QC-filtered long reads
 * - `fna?`: Assembly file (passed through for assembly-based samples)
 */
nextflow.preview.types = true

include { QC as QC_MODULE } from '../../../modules/bactopia/qc/main'
include { filterWithData  } from 'plugin/nf-bactopia'

workflow QC {
    take:
    samples: Channel<Record>
    adapters: Path?
    phix: Path?

    main:
    ch_qc = QC_MODULE(samples, adapters, phix)

    emit:
    // Downstream inputs
    reads = filterWithData(ch_qc, ['r1', 'r2', 'se', 'lr', 'fna'])

    // Published outputs
    sample_outputs = ch_qc
    run_outputs = channel.empty()
}
