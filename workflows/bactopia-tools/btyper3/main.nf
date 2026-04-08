#!/usr/bin/env nextflow
/**
 * Taxonomic classification of Bacillus cereus group isolates.
 *
 * This Bactopia Tool uses [BTyper3](https://github.com/lmc297/BTyper3) to classify
 * Bacillus cereus group isolates from genome assemblies.
 *
 * @status stable
 * @keywords bacillus cereus, taxonomy, classification, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,classification
 * @citation btyper3
 *
 * @subworkflows bactopiatool_init, btyper3
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *_final_results.txt   Final tab-delimited file of BTyper3 results
 * @publish results/*             Directory of detailed analysis results
 *
 * @section Merged Results
 * @publish btyper3.tsv           Merged TSV file containing BTyper3 results from all samples
 *
 * @section Execution Logs
 * @publish logs/**               Tool execution logs
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BTYPER3             } from '../../../subworkflows/btyper3/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_btyper3 = BTYPER3(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_btyper3.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_btyper3.sample_outputs)
    // Run-level
    run_outputs = ch_btyper3.run_outputs
    run_nf_logs = collectNextflowLogs(ch_btyper3.run_outputs)
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
