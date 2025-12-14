#!/usr/bin/env nextflow
/**
 * Removal of human and contaminant sequences from metagenomic reads.
 *
 * This Bactopia Tool removes human and other contaminant sequences from metagenomic reads using
 * either [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) or a Kraken2-based approach
 * (k2scrubber) with the HPRC human database. The tool provides flexible contamination removal
 * with detailed reporting of read classification and filtering statistics. It processes paired-end
 * or single-end reads, producing cleaned FASTQ files with human sequences removed and comprehensive
 * reports documenting the decontamination process.
 *
 * @status stable
 * @keywords metagenomics, decontamination, human removal, read filtering, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation,conditional
 * @citation kraken2, srahumanscrubber
 *
 * @subworkflows bactopiatool_init, scrubber
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input use_srascrubber
 * Boolean flag to choose between SRA Human Scrubber (true) or k2scrubber (false) for decontamination
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
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    use_srascrubber : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SCRUBBER          } from '../../../subworkflows/scrubber/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    SCRUBBER(
        BACTOPIATOOL_INIT.out.reads,
        params.use_srascrubber
    )

    // Collect outputs
    ch_results = ch_results.mix(SCRUBBER.out.results)
    ch_logs = ch_logs.mix(SCRUBBER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SCRUBBER.out.nf_logs)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)

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
