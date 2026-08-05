#!/usr/bin/env nextflow
/**
 * Identify and type Stx operons from assembled genomic sequences
 *
 * This Bactopia Tool uses [StxTyper](https://github.com/ncbi/stxtyper) to identify and type stx operons from assembled genomic sequences.
 *
 * @status stable
 * @keywords stx, shiga toxin, typing, stec, virulence, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, stxtyper
 *
 * @subworkflows utils_bactopia-tools, stxtyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                        Tab-delimited Stx operon typing results
 *
 * @section Merged Results
 * @publish stxtyper.tsv              Merged TSV file containing stxtyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/stxtyper/*               Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { STXTYPER            } from '../../../subworkflows/stxtyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_stxtyper = STXTYPER(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_stxtyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_stxtyper.sample_outputs)
    // Run-level
    run_outputs = ch_stxtyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_stxtyper.run_outputs)
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
