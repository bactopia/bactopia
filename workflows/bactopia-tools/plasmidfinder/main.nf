#!/usr/bin/env nextflow
/**
 * Bactopia Tool: Plasmidfinder.
 *
 * Plasmid identification from assemblies
 * The `plasmidfinder` module identifies plasmids in total or partial sequenced isolates of bacteria.
 *
 * @status stable
 * @keywords plasmid, identification, replicon, typing, assembly
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation plasmidfinder
 *
 * @subworkflows utils_bactopia-tools, plasmidfinder
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @section Per-Sample Results
 * @publish *    Analysis results
 *
 * @section Merged Results
 * @publish merged-*    Aggregated results from all samples
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml Software version information
  */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { PLASMIDFINDER       } from '../../../subworkflows/plasmidfinder/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_plasmidfinder = PLASMIDFINDER(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_plasmidfinder.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_plasmidfinder.sample_outputs)
    // Run-level
    run_outputs = ch_plasmidfinder.run_outputs
    run_nf_logs = collectNextflowLogs(ch_plasmidfinder.run_outputs)
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
