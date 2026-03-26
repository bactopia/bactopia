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
 * @subworkflows bactopiatool_init, mashdist
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
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    mash_sketch : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MASHDIST          } from '../../../subworkflows/mashdist/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MASHDIST(
        BACTOPIATOOL_INIT.out.assembly,
        params.mash_sketch
    )
    ch_sample_nf_logs = MASHDIST.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MASHDIST.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = MASHDIST.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MASHDIST.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
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
