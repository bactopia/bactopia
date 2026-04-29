#!/usr/bin/env nextflow
/**
 * Run all available MLST schemes for a species against an assembly
 *
 * This Bactopia Tool uses [GigaTyper](https://github.com/rpetit3/gigatyper) to run all available mlst schemes for a species against an assembly.
 *
 * @status stable
 * @keywords mlst, typing, multi-scheme, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, gigatyper
 *
 * @subworkflows utils_bactopia-tools, gigatyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                     MLST results across all schemes
 *
 * @section Merged Results
 * @publish gigatyper.tsv             Merged TSV file containing gigatyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/gigatyper/*          Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                 Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml              Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { GIGATYPER           } from '../../../subworkflows/gigatyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_gigatyper = GIGATYPER(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_gigatyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_gigatyper.sample_outputs)
    // Run-level
    run_outputs = ch_gigatyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_gigatyper.run_outputs)
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
