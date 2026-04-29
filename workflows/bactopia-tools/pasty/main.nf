#!/usr/bin/env nextflow
/**
 * In silico serogrouping of Pseudomonas aeruginosa isolates.
 *
 * This Bactopia Tool uses [pasty](https://github.com/rpetit3/pasty) for
 * serogrouping of _Pseudomonas aeruginosa_ isolates from genome assemblies.
 *
 * @status stable
 * @keywords fasta, serogrouping, Pseudomonas aeruginosa, typing, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation pasty
 *
 * @subworkflows utils_bactopia-tools, pasty
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Serogrouping results for each sample
 *
 * @section Merged Results
 * @publish pasty.tsv        Merged TSV file containing Pasty results from all samples
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
include { PASTY               } from '../../../subworkflows/pasty/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_pasty = PASTY(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_pasty.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_pasty.sample_outputs)
    // Run-level
    run_outputs = ch_pasty.run_outputs
    run_nf_logs = collectNextflowLogs(ch_pasty.run_outputs)
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
