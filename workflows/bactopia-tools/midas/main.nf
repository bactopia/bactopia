#!/usr/bin/env nextflow
/**
 * Estimate species abundances from metagenomic samples.
 *
 * This Bactopia Tool uses [MIDAS](https://github.com/snayfach/MIDAS) to estimate
 * bacterial species abundances in metagenomic samples. MIDAS uses a database
 * with more than 30,000 reference genomes for accurate species profiling.
 *
 * @status stable
 * @keywords metagenomics, species abundance, profiling, midas, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,profiling,database-dependent
 * @citation midas
 *
 * @subworkflows bactopiatool_init, midas
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input midas_db
 * Path to MIDAS database
 *
 * @input download_midas
 * Download MIDAS database if not found locally
 *
 * @section Species Abundance
 * @publish *.tsv              Species abundance profiles
 * @publish *-species.tsv       Species-level abundance
 * @publish *-genes.tsv         Gene-level abundance
 *
 * @section Merged Results
 * @publish midas.tsv           Merged TSV file containing MIDAS results from all samples
 *
 * @section Execution Logs
 * @publish logs/**             Tool execution logs
 * @publish logs/nf-*           Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    midas_db       : Path
    download_midas : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MIDAS             } from '../../../subworkflows/midas/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MIDAS(
        BACTOPIATOOL_INIT.out.reads,
        params.midas_db,
        params.download_midas,
        params.midas_save_as_tarball
    )
    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = MIDAS.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MIDAS.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = MIDAS.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = MIDAS.out.run_outputs
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
