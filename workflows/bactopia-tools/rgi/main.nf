#!/usr/bin/env nextflow
/**
 * Prediction of antibiotic resistance genes using RGI.
 *
 * This Bactopia Tool uses [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) to identify
 * and characterize antibiotic resistance genes in bacterial assemblies. RGI integrates with the
 * Comprehensive Antibiotic Resistance Database (CARD) to provide high-confidence predictions
 * of resistance determinants, including perfect and strict hits to known resistance genes,
 * as well as loose hits for novel variants. The tool generates detailed reports in both
 * JSON and TSV formats, along with heatmap visualizations for comparative analysis
 * across multiple samples.
 *
 * @status stable
 * @keywords bacteria, antibiotic resistance, card, resistance genes, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation,visualization
 * @citation csvtk, rgi
 *
 * @subworkflows bactopiatool_init, rgi
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.json                    JSON report containing detailed RGI results for each sample
 * @publish *.txt                     Tab-delimited report of RGI results for each sample
 *
 * @section Merged Results
 * @publish rgi.tsv                   Merged TSV file containing RGI results from all samples
 * @publish rgi-2.{csv,png}           Heatmap representations of resistance genes across all samples
 *
 * @section Execution Logs
 * @publish logs/**                    Tool execution logs including RGI output
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { RGI               } from '../../../subworkflows/rgi/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    RGI(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(RGI.out.results)
    ch_logs = ch_logs.mix(RGI.out.logs)
    ch_nf_logs = ch_nf_logs.mix(RGI.out.nf_logs)
    ch_versions = ch_versions.mix(RGI.out.versions)

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
