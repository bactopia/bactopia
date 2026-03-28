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
 * @subworkflows bactopiatool_init, pneumocat
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
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PNEUMOCAT         } from '../../../subworkflows/pneumocat/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PNEUMOCAT(BACTOPIATOOL_INIT.out.reads)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = PNEUMOCAT.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = PNEUMOCAT.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = PNEUMOCAT.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = PNEUMOCAT.out.run_outputs
    run_nf_logs = ch_run_nf_logs
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
