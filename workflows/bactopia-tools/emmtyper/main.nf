#!/usr/bin/env nextflow
/**
 * emm-typing of Streptococcus pyogenes assemblies.
 *
 * This Bactopia Tool uses [emmtyper](https://github.com/MDU-PHL/emmtyper) for
 * emm-typing of _Streptococcus pyogenes_ using a de novo or complete assembly.
 *
 * @status stable
 * @keywords streptococcus pyogenes, emm typing, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation emmtyper
 *
 * @subworkflows bactopiatool_init, emmtyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input emmtyper_blastdb
 * Path to emmtyper BLAST database
 *
 * @section Per-Sample Results
 * @publish *.tsv           Tab-delimited file with emmtyper result
 *
 * @section Merged Results
 * @publish emmtyper.tsv     Merged TSV file containing emmtyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    emmtyper_blastdb : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { EMMTYPER          } from '../../../subworkflows/emmtyper/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    EMMTYPER(
        BACTOPIATOOL_INIT.out.assembly,
        params.emmtyper_blastdb
    )
    ch_sample_nf_logs = EMMTYPER.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = EMMTYPER.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = EMMTYPER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = EMMTYPER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}
output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
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

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
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
