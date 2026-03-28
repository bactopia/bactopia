#!/usr/bin/env nextflow
/**
 * Taxonomic profiling by abundance-corrected MinHash.
 *
 * This Bactopia Tool uses [Sylph](https://github.com/bluenote-1577/Sylph) to perform
 * taxonomic profiling of metagenomic samples using abundance-corrected MinHash sketches
 * for accurate species-level quantification.
 *
 * @status stable
 * @keywords taxonomic profiling, metagenomics, minhash, abundance, sylph, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,profiling
 * @citation sylph
 *
 * @subworkflows bactopiatool_init, sylph
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input sylph_db
 * Path to Sylph database for taxonomic classification
 *
 * @section Per-Sample Results
 * @publish *.profile.txt       Species abundance profile
 * @publish *.krona.html        Interactive Krona plot visualization
 *
 * @section Merged Results
 * @publish sylph.tsv          Merged TSV file containing Sylph profiles from all samples
 *
 * @section Execution Logs
 * @publish logs/**            Tool execution logs
 * @publish logs/nf-*          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml       Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    sylph_db : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SYLPH             } from '../../../subworkflows/sylph/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SYLPH(BACTOPIATOOL_INIT.out.reads, params.sylph_db)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = SYLPH.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = SYLPH.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = SYLPH.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = SYLPH.out.run_outputs
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
