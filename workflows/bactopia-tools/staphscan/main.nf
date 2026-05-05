#!/usr/bin/env nextflow
/**
 * Genome-based surveillance analysis of Staphylococcus aureus.
 *
 * This Bactopia Tool uses [StaphSCAN](https://github.com/riccabolla/StaphSCAN) to perform
 * genome-based surveillance of _Staphylococcus aureus_ for epidemiological typing and
 * resistance profiling.
 *
 * @status stable
 * @keywords staphylococcus aureus, surveillance, mlst, spa typing, sccmec, amr, virulence, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, staphscan
 *
 * @subworkflows utils_bactopia-tools, staphscan
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input staphscan_db_mlst
 * Path or tarball to custom MLST database (optional)
 *
 * @section Per-Sample Results
 * @publish *.tsv                        Per-sample surveillance summary with MLST, spa type, SCCmec, capsule, AGR, resistance, biofilm, and virulence results
 *
 * @section Merged Results
 * @publish staphscan.tsv             Merged TSV file containing staphscan results from all samples
 *
 * @section Execution Logs
 * @publish logs/staphscan/*              Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    staphscan_db_mlst : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { STAPHSCAN           } from '../../../subworkflows/staphscan/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_staphscan = STAPHSCAN(ch_bactopiatool.assembly, params.staphscan_db_mlst)

    publish:
    // Per-sample
    sample_outputs = ch_staphscan.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_staphscan.sample_outputs)
    // Run-level
    run_outputs = ch_staphscan.run_outputs
    run_nf_logs = collectNextflowLogs(ch_staphscan.run_outputs)
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
