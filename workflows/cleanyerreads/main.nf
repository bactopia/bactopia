#!/usr/bin/env nextflow
/**
 * Quality control and optional host read removal from raw sequencing reads.
 *
 * This workflow performs comprehensive read quality control including trimming,
 * adapter removal, quality filtering, and optionally removes host contamination
 * using [Kraken2](https://github.com/DerrickWood/kraken2) or SRA Scrubber.
 * It processes raw sequencing reads to produce high-quality clean reads ready
 * for downstream analysis.
 *
 * @status stable
 * @keywords reads, quality control, trimming, filtering, host removal, preprocessing
 * @tags complexity:moderate input-type:parameter output-type:multiple features:aggregation,conditional-logic,database-dependent
 *
 * @subworkflows bactopia_init, gather, qc, scrubber
 *
 * @input rundir
 * Directory containing raw sequencing reads
 *
 * @input use_k2scrubber
 * Remove host reads using Kraken2 scrubber
 *
 * @input use_srascrubber
 * Remove host reads using SRA scrubber
 *
 * @input adapters
 * Path to adapter sequences file for removal
 *
 * @input phix
 * Path to PhiX sequences for contamination removal
 *
 * @section Quality Control Reports
 * @publish supplemental/*_fastqc.*          FastQC reports for raw and cleaned reads
 * @publish supplemental/*-NanoPlot.*       NanoPlot reports for Nanopore reads
 * @publish supplemental/*.fastp.*          Fastp quality reports (when applicable)
 * @publish supplemental/*_original.json    Quality metrics for original reads
 * @publish supplemental/*_final.json       Quality metrics for final reads
 *
 * @section Cleaned Reads
 * @publish *.fastq.gz                       Quality controlled and trimmed reads
 * @publish *.unclassified.fastq.gz          Host-decontaminated reads (if scrubber enabled)
 *
 * @section Host Decontamination
 * @note Only created when host removal is enabled
 * @publish *.kraken2.report.txt             Kraken2 classification report
 * @publish *.scrub.report.tsv               Human contamination report
 *
 * @section Execution Logs
 * @publish logs/**                          Tool execution logs
 * @publish logs/nf-*                        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                     Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    use_k2scrubber  : Boolean
    use_srascrubber : Boolean
    adapters        : Path?
    phix            : Path?
}

// Core
include { BACTOPIA_INIT   } from '../../subworkflows/utils/bactopia/main'
include { GATHER          } from '../../subworkflows/bactopia/gather/main'
include { QC              } from '../../subworkflows/bactopia/qc/main'

// Scrubber
include { SCRUBBER        } from '../../subworkflows/scrubber/main'

include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIA_INIT()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)

    // Collect outputs
    ch_sample_outputs = GATHER.out.sample_outputs
    ch_run_outputs = GATHER.out.run_outputs

    if (params.use_k2scrubber || params.use_srascrubber) {
        // Remove host reads
        SCRUBBER(GATHER.out.reads, params.use_srascrubber, params.nohuman_db ? file(params.nohuman_db) : file("NO_DB"), params.download_nohuman, params.nohuman_save_as_tarball)
        ch_sample_outputs = ch_sample_outputs.mix(SCRUBBER.out.sample_outputs)
        ch_run_outputs = ch_run_outputs.mix(SCRUBBER.out.run_outputs)

        // Clean up scrubbed reads
        QC(SCRUBBER.out.scrubbed_extra, params.adapters, params.phix)
    } else {
        // Clean up raw reads
        QC(GATHER.out.reads, params.adapters, params.phix)
    }
    ch_sample_outputs = ch_sample_outputs.mix(QC.out.sample_outputs)

    ch_sample_nf_logs = collectNextflowLogs(ch_sample_outputs)
    ch_run_nf_logs = collectNextflowLogs(ch_run_outputs)

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ch_sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ch_run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
