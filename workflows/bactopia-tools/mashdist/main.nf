#!/usr/bin/env nextflow
/**
 * Calculate Mash distances between sequences and reference genomes.
 *
 * This Bactopia Tool uses [Mash](https://github.com/marbl/Mash) to determine the Mash
 * distance from samples to reference genome sketches for rapid genomic comparison.
 *
 * @status stable
 * @keywords mash, distance, similarity, comparative genomics, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,comparative
 * @citation mash
 *
 * @subworkflows utils_bactopia-tools, mashdist
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input mash_sketch
 * Path to Mash sketch file of reference genome(s)
 *
 * @section Per-Sample Results
 * @publish *.txt            Mash distance results for each sample
 *
 * @section Merged Results
 * @publish mashdist.tsv      Merged TSV file containing Mash distances from all samples
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
    mash_sketch : Path
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { MASHDIST            } from '../../../subworkflows/mashdist/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_mashdist = MASHDIST(ch_bactopiatool.assembly, params.mash_sketch)

    publish:
    // Per-sample
    sample_outputs = ch_mashdist.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_mashdist.sample_outputs)
    // Run-level
    run_outputs = ch_mashdist.run_outputs
    run_nf_logs = collectNextflowLogs(ch_mashdist.run_outputs)
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
