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
 * @subworkflows utils_bactopia-tools, rgi
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
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { RGI                 } from '../../../subworkflows/rgi/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_rgi = RGI(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_rgi.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_rgi.sample_outputs)
    // Run-level
    run_outputs = ch_rgi.run_outputs
    run_nf_logs = collectNextflowLogs(ch_rgi.run_outputs)
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
