#!/usr/bin/env nextflow
/**
 * Bactopia Tool: Plasmidfinder.
 *
 * Plasmid identification from assemblies
 * The `plasmidfinder` module identifies plasmids in total or partial sequenced isolates of bacteria.
 *
 * @status stable
 *
 * @subworkflows bactopiatool_init, plasmidfinder
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
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PLASMIDFINDER     } from '../../../subworkflows/plasmidfinder/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PLASMIDFINDER(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = PLASMIDFINDER.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = PLASMIDFINDER.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = PLASMIDFINDER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = PLASMIDFINDER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}
output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
