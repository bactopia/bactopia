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
 * @subworkflows utils_bactopia-tools, midas
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
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    midas_db       : Path
    download_midas : Boolean
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { MIDAS               } from '../../../subworkflows/midas/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_midas = MIDAS(
        ch_bactopiatool.reads,
        params.midas_db,
        params.download_midas,
        params.midas_save_as_tarball
    )

    publish:
    // Per-sample
    sample_outputs = ch_midas.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_midas.sample_outputs)
    // Run-level
    run_outputs = ch_midas.run_outputs
    run_nf_logs = collectNextflowLogs(ch_midas.run_outputs)
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
