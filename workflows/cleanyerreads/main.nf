#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Clean Yer Reads.
 *
 * Quality control and optional host read removal from raw sequencing reads.
 * This workflow performs read trimming, adapter removal, quality filtering,
 * and optionally removes host contamination using Kraken2 or SRA Scrubber.
 *
 * @status stable
 * @keywords reads, quality control, trimming, filtering, host removal
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
 * @section Quality Reports
 * @publish fastqc/*   FastQC quality control reports
 * @publish multiqc/*   MultiQC aggregated quality reports
 *
 * @section Cleaned Reads
 * @publish *.fastq.gz   Quality controlled and trimmed reads
 * @publish *.scrubbed.fastq.gz   Host-decontaminated reads (if scrubber used)
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution logs
 *
 * @section Versions
 * @publish versions.yml Software version information
   */

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

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIA_INIT()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)
    ch_results = ch_results.mix(GATHER.out.results)
    ch_logs = ch_logs.mix(GATHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
    ch_versions = ch_versions.mix(GATHER.out.versions)

    if (params.use_k2scrubber || params.use_srascrubber) {
        // Remove host reads
        SCRUBBER(GATHER.out.fastq_only, params.use_srascrubber)
        ch_results = ch_results.mix(SCRUBBER.out.results)
        ch_logs = ch_logs.mix(SCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SCRUBBER.out.nf_logs)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Clean up scrubbed reads
        QC(
            SCRUBBER.out.scrubbed_extra,
            params.adapters,
            params.phix
        )
    } else {
        // Clean up raw reads
        QC(
            GATHER.out.raw_fastq,
            params.adapters,
            params.phix
        )
    }
    ch_results = ch_results.mix(QC.out.results)
    ch_logs = ch_logs.mix(QC.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QC.out.nf_logs)
    ch_versions = ch_versions.mix(QC.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
