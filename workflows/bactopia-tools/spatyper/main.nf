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
 * @subworkflows bactopiatool_init, spatyper
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
    rundir   : String

    // Tool-specific parameters
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SPATYPER          } from '../../../subworkflows/spatyper/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SPATYPER(
        BACTOPIATOOL_INIT.out.assembly,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )
    ch_sample_nf_logs = SPATYPER.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = SPATYPER.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = SPATYPER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = SPATYPER.out.run_outputs
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
