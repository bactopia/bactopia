#!/usr/bin/env nextflow
/**
 * In silico serogrouping of Pseudomonas aeruginosa isolates.
 *
 * This Bactopia Tool uses [pasty](https://github.com/rpetit3/pasty) for
 * serogrouping of _Pseudomonas aeruginosa_ isolates from genome assemblies.
 *
 * @status stable
 * @keywords fasta, serogrouping, Pseudomonas aeruginosa, typing, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation pasty
 *
 * @subworkflows bactopiatool_init, pasty
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Serogrouping results for each sample
 *
 * @section Merged Results
 * @publish pasty.tsv        Merged TSV file containing Pasty results from all samples
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
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PASTY             } from '../../../subworkflows/pasty/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PASTY(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = PASTY.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = PASTY.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = PASTY.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = PASTY.out.run_outputs
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
