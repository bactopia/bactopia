#!/usr/bin/env nextflow
/**
 * Sequence Based Typing (SBT) of Legionella pneumophila.
 *
 * This Bactopia Tool uses [legsta](https://github.com/tseemann/legsta) for
 * _in silico_ _Legionella pneumophila_ Sequence Based Typing (SBT).
 *
 * @status stable
 * @keywords legionella pneumophila, sbt, typing, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation legsta
 *
 * @subworkflows bactopiatool_init, legsta
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            SBT typing results
 *
 * @section Merged Results
 * @publish legsta.tsv        Merged TSV file containing legsta results from all samples
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

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { LEGSTA              } from '../../../subworkflows/legsta/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_legsta = LEGSTA(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_legsta.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_legsta.sample_outputs)
    // Run-level
    run_outputs = ch_legsta.run_outputs
    run_nf_logs = collectNextflowLogs(ch_legsta.run_outputs)
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
