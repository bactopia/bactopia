#!/usr/bin/env nextflow
/**
 * Serovar prediction of Salmonella enterica from assemblies.
 *
 * This Bactopia Tool uses [Salmonella In Silico Typing Resource](https://github.com/phac-nml/sistr_cmd),
 * or SISTR, for serovar prediction of Salmonella enterica assemblies using cgMLST typing
 * and molecular serovar prediction.
 *
 * @status stable
 * @keywords salmonella, serovar, cgmlst, typing, sistr, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing,cgmlst
 * @citation sistr
 *
 * @subworkflows bactopiatool_init, sistr
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.csv            SISTR analysis results in CSV format
 *
 * @section Merged Results
 * @publish sistr.tsv         Merged TSV file containing SISTR results from all samples
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
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SISTR             } from '../../../subworkflows/sistr/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SISTR(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = SISTR.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = SISTR.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = SISTR.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = SISTR.out.run_outputs
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
