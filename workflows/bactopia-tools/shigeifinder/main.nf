#!/usr/bin/env nextflow
/**
 * In silico serotype prediction for Shigella and Enteroinvasive E. coli (EIEC).
 *
 * This Bactopia Tool uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to conduct rapid
 * serotype prediction from genome assemblies. It provides species identification and predicts
 * the traditional O and H antigens for Shigella and EIEC isolates, enabling epidemiological
 * tracking and surveillance without requiring traditional serological methods.
 *
 * @status stable
 * @keywords Shigella, EIEC, serotyping, prediction, epidemiology
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, shigeifinder
 *
 * @subworkflows bactopiatool_init, shigeifinder
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                      Tab-delimited file containing predicted Shigella or EIEC serotype
 *
 * @section Merged Results
 * @publish shigeifinder.tsv            Merged TSV file containing serotype predictions from all samples
 *
 * @section Execution Logs
 * @publish logs/shigeifinder/*         Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                   Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SHIGEIFINDER      } from '../../../subworkflows/shigeifinder/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SHIGEIFINDER(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = SHIGEIFINDER.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = SHIGEIFINDER.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = SHIGEIFINDER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = SHIGEIFINDER.out.run_outputs
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
