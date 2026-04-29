#!/usr/bin/env nextflow
/**
 * Removal of human and contaminant sequences from metagenomic reads.
 *
 * This Bactopia Tool removes human and other contaminant sequences from metagenomic reads using
 * either [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) or
 * [nohuman](https://github.com/mbhall88/nohuman) with the HPRC human database. The tool provides flexible contamination removal
 * with detailed reporting of read classification and filtering statistics. It processes paired-end
 * or single-end reads, producing cleaned FASTQ files with human sequences removed and comprehensive
 * reports documenting the decontamination process.
 *
 * @status stable
 * @keywords metagenomics, decontamination, human removal, read filtering, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation,conditional
 * @citation kraken2, srahumanscrubber
 *
 * @subworkflows utils_bactopia-tools, scrubber
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input use_srascrubber
 * Boolean flag to choose between SRA Human Scrubber (true) or nohuman (false) for decontamination
 *
 * @section Per-Sample Results
 * @publish *.scrubbed.fastq.gz          Cleaned reads after human sequence removal
 * @publish *.scrub.report.tsv           Report of read classification and removal statistics
 *
 * @section Merged Results
 * @publish scrubber.tsv                 Merged TSV file containing scrubber reports from all samples
 *
 * @section Execution Logs
 * @publish logs/**                      Tool execution logs including classification output
 * @publish logs/nf-*                    Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                 Software version information
 */
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    use_srascrubber         : Boolean
    nohuman_db              : Path?
    download_nohuman        : Boolean
    nohuman_save_as_tarball : Boolean
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SCRUBBER            } from '../../../subworkflows/scrubber/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_scrubber = SCRUBBER(
        ch_bactopiatool.reads,
        params.use_srascrubber,
        params.nohuman_db,
        params.download_nohuman,
        params.nohuman_save_as_tarball
    )

    publish:
    // Per-sample
    sample_outputs = ch_scrubber.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_scrubber.sample_outputs)
    // Run-level
    run_outputs = ch_scrubber.run_outputs
    run_nf_logs = collectNextflowLogs(ch_scrubber.run_outputs)
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
