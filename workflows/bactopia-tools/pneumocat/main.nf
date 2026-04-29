#!/usr/bin/env nextflow
/**
 * Capsular type assignment to Streptococcus pneumoniae from sequence reads.
 *
 * This Bactopia Tool uses [PneumoCaT](https://github.com/ukhsa-collaboration/PneumoCaT) to assign capsular
 * type to _Streptococcus pneumoniae_ from sequence reads for epidemiological typing.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, capsular typing, pneumocat, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation pneumocat
 *
 * @subworkflows utils_bactopia-tools, pneumocat
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Capsular type assignment results
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { PNEUMOCAT           } from '../../../subworkflows/pneumocat/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_pneumocat = PNEUMOCAT(ch_bactopiatool.reads)

    publish:
    // Per-sample
    sample_outputs = ch_pneumocat.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_pneumocat.sample_outputs)
    // Run-level
    run_outputs = ch_pneumocat.run_outputs
    run_nf_logs = collectNextflowLogs(ch_pneumocat.run_outputs)
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
