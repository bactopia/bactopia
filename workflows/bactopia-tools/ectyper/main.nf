#!/usr/bin/env nextflow
/**
 * In silico prediction of Escherichia coli serotype.
 *
 * This Bactopia Tool uses [ECTyper](https://github.com/phac-nml/ecoli_serotyping) to conduct
 * _in silico_ prediction of serotype for _Escherichia coli_ genomes. It uses the genome assemblies
 * to provide basic species identification and the predicted _E. coli_ serotype (e.g. O174:H21).
 *
 * @status stable
 * @keywords escherichia coli, serotyping, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation ectyper
 *
 * @subworkflows bactopiatool_init, ectyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                  Tab-delimited file with ECTyper result
 * @publish blast_output_alleles.txt Allele report generated from BLAST results
 *
 * @section Merged Results
 * @publish ectyper.tsv            Merged TSV file containing ECTyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                 Tool execution logs
 * @publish logs/nf-*               Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ECTYPER           } from '../../../subworkflows/ectyper/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    ECTYPER(BACTOPIATOOL_INIT.out.assembly)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ECTYPER.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = ECTYPER.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ECTYPER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ECTYPER.out.run_outputs
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
