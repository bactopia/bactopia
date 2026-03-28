#!/usr/bin/env nextflow
/**
 * Identification, classification, and annotation of translated gene matches.
 *
 * This Bactopia Tool uses [GAMMA](https://github.com/rastanton/GAMMA) to identify, classify, and annotate
 * translated gene matches from assemblies using a comprehensive protein database.
 *
 * @status stable
 * @keywords gene annotation, protein classification, translation, gamma, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,annotation
 * @citation gamma
 *
 * @subworkflows bactopiatool_init, gamma
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input gamma_db
 * Path to GAMMA database
 *
 * @section Per-Sample Results
 * @publish *.tsv            Tab-delimited file with gene classification results
 *
 * @section Merged Results
 * @publish gamma.tsv         Merged TSV file containing GAMMA results from all samples
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
    gamma_db : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { GAMMA             } from '../../../subworkflows/gamma/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    GAMMA(
        BACTOPIATOOL_INIT.out.assembly,
        params.gamma_db
    )

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = GAMMA.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = GAMMA.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = GAMMA.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = GAMMA.out.run_outputs
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
