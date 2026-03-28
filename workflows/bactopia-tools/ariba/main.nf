#!/usr/bin/env nextflow
/**
 * Gene identification through local assemblies.
 *
 * This Bactopia Tool uses [ARIBA](https://github.com/sanger-pathogens/ariba) to rapidly
 * identify genes in a database by creating local assemblies from short-read data.
 * ARIBA performs reference-based assembly and variant calling for gene detection.
 *
 * @status stable
 * @keywords fastq, assembly, resistance, virulence, gene detection, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,database-dependent,gene-identification
 * @citation ariba
 *
 * @subworkflows bactopiatool_init, ariba
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input ariba_db
 * Name of the ARIBA database to use for gene identification
 *
 * @section Per-Sample Results
 * @publish *-report.tsv           Gene detection report for each sample
 * @publish *-summary.csv          Summary of gene detection results
 * @publish assembled_genes.fa.gz  Assembled genes in compressed FASTA format
 * @publish assembled_seqs.fa.gz   Assembled sequences matching references
 * @publish assemblies.fa.gz       Raw local assemblies
 * @publish debug.report.tsv       Detailed report including synonymous mutations
 * @publish log.clusters.gz        Analysis log file
 * @publish version_info.txt       Version information for ARIBA and dependencies
 *
 * @section Merged Results
 * @publish ariba-report.tsv       Merged gene detection reports from all samples
 * @publish ariba-summary.csv      Merged summary reports from all samples
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml           Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    ariba_db : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ARIBA             } from '../../../subworkflows/ariba/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    ARIBA(BACTOPIATOOL_INIT.out.reads, params.ariba_db)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ARIBA.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = ARIBA.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ARIBA.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ARIBA.out.run_outputs
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
