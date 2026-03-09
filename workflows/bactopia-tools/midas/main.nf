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
    midas_db : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MIDAS             } from '../../../subworkflows/midas/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MIDAS(
        BACTOPIATOOL_INIT.out.reads,
        params.midas_db
    )
    ch_sample_nf_logs = MIDAS.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MIDAS.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = MIDAS.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MIDAS.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
