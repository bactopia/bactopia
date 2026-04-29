#!/usr/bin/env nextflow
/**
 * Serotype prediction of Haemophilus parasuis assemblies.
 *
 * This Bactopia Tool uses [HpsuisSero](https://github.com/jimmyliu1326/HpsuisSero) to predict
 * the serotype of _Haemophilus parasuis_ assemblies.
 *
 * @status stable
 * @keywords haemophilus parasuis, serotyping, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation hpsuissero
 *
 * @subworkflows utils_bactopia-tools, hpsuissero
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.summary       Summary of serotype prediction
 *
 * @section Merged Results
 * @publish hpsuissero.tsv   Merged TSV file containing HpsuisSero results from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml       Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { HPSUISSERO          } from '../../../subworkflows/hpsuissero/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_hpsuissero = HPSUISSERO(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_hpsuissero.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_hpsuissero.sample_outputs)
    // Run-level
    run_outputs = ch_hpsuissero.run_outputs
    run_nf_logs = collectNextflowLogs(ch_hpsuissero.run_outputs)
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
