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
 * @subworkflows utils_bactopia-tools, shigeifinder
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

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SHIGEIFINDER        } from '../../../subworkflows/shigeifinder/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_shigeifinder = SHIGEIFINDER(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_shigeifinder.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_shigeifinder.sample_outputs)
    // Run-level
    run_outputs = ch_shigeifinder.run_outputs
    run_nf_logs = collectNextflowLogs(ch_shigeifinder.run_outputs)
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
