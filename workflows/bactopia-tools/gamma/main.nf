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
 * @subworkflows utils_bactopia-tools, gamma
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
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    gamma_db : Path
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { GAMMA               } from '../../../subworkflows/gamma/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_gamma = GAMMA(ch_bactopiatool.assembly, params.gamma_db)

    publish:
    // Per-sample
    sample_outputs = ch_gamma.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_gamma.sample_outputs)
    // Run-level
    run_outputs = ch_gamma.run_outputs
    run_nf_logs = collectNextflowLogs(ch_gamma.run_outputs)
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
