#!/usr/bin/env nextflow
/**
 * spa typing of Staphylococcus aureus assemblies.
 *
 * This Bactopia Tool uses [spaTyper](https://github.com/HCGB-IGTP/spaTyper) to assign
 * spa types to _Staphylococcus aureus_ assemblies for epidemiological typing.
 *
 * @status stable
 * @keywords staphylococcus aureus, spa typing, epidemiology, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation spatyper
 *
 * @subworkflows utils_bactopia-tools, spatyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input spatyper_repeats
 * Path to spaTyper repeats database
 *
 * @input spatyper_repeat_order
 * Path to spaTyper repeat order file
 *
 * @section Per-Sample Results
 * @publish *.txt            spa typing results for each sample
 *
 * @section Merged Results
 * @publish spatyper.tsv      Merged TSV file containing spaTyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SPATYPER            } from '../../../subworkflows/spatyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_spatyper = SPATYPER(
        ch_bactopiatool.assembly,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )

    publish:
    // Per-sample
    sample_outputs = ch_spatyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_spatyper.sample_outputs)
    // Run-level
    run_outputs = ch_spatyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_spatyper.run_outputs)
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
